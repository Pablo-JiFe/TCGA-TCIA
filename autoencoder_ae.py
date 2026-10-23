import pandas as pd
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import TensorDataset, DataLoader
from torch_snippets import *
import random

random.seed(37)

device = "cuda" if torch.cuda.is_available() else "cpu"

counts = pd.read_parquet("./data/counts_disease5k_tcga.parquet")
counts.rename(counts["gene_name"], inplace=True)
counts = counts.drop(columns=["gene_name"])
counts = counts.T

counts_train = counts.sample(n=700)
counts_test = counts.loc[~ counts.index.isin(counts_train.index)]


def normalize(x):
    mean = x.mean()
    sd = x.std()
    return (x - mean) / (sd + 1e-8)

counts_train_norm = normalize(counts_train)
counts_test_norm = normalize(counts_test)
    
counts_ten_train = torch.from_numpy(counts_train_norm.values).float()
counts_ten_test = torch.from_numpy(counts_test_norm.values).float()

train_ds = TensorDataset(counts_ten_train)
test_ds = TensorDataset(counts_ten_test)

train_loader = DataLoader(train_ds, batch_size=64, shuffle=True)
test_loader = DataLoader(test_ds, batch_size=64, shuffle=False)


class AutoEncoder(nn.Module):
    def __init__(self, input_dim = 5000, latent_dim = 5, hidden_deep = 750, hidden_sup = 3000):
        super().__init__()
        self.encoder = nn.Sequential(
            nn.Linear(input_dim, hidden_sup),
            nn.LayerNorm(hidden_sup),
            nn.GELU(),
            nn.Linear(hidden_sup, hidden_deep),
            nn.LayerNorm(hidden_deep),
            nn.GELU(),
            nn.Linear(hidden_deep, latent_dim)
        )
        self.decoder = nn.Sequential(
            nn.Linear(latent_dim, hidden_deep),
            nn.LayerNorm(hidden_deep),
            nn.GELU(),
            nn.Linear(hidden_deep, hidden_sup),
            nn.LayerNorm(hidden_sup),
            nn.GELU(),
            nn.Linear(hidden_sup, input_dim)
        )
        
    def forward(self, x):
        to_latent = self.encoder(x)
        decode = self.decoder(to_latent)
        return decode

model = AutoEncoder().to(device)
optimizer = optim.Adam(model.parameters(), lr=1e-3, weight_decay=1e-3)
criterion = nn.HuberLoss()
epochs = 30


def train_batch(batch, model, criterion, optimizer):
    model.train()
    optimizer.zero_grad()
    input = batch[0].to(device)
    output = model(input)
    loss = criterion(output, input)
    loss.backward()
    optimizer.step()
    return loss.item()


@torch.no_grad()
def validate_batch(batch, model, criterion):
    model.eval()
    input = batch[0].to(device)
    output = model(input)
    loss = criterion(output, input)
    return loss.item()


log = Report(epochs)

for epoch in range(epochs):

    N = len(train_loader)
    train_loss = 0
    for batch in train_loader:
        loss = train_batch(batch, model, criterion, optimizer)
        train_loss += loss
    train_loss /= N

    val_loss = 0
    for batch in test_loader:
        val_loss += validate_batch(batch, model, criterion)
    val_loss /= len(test_loader)

    log.record(pos=epoch + 1, trn_loss=train_loss, val_loss=val_loss, end='\r')
    log.report_avgs(epoch + 1)

log.plot_epochs(log=True)

model.eval()

with torch.no_grad():
    X_train = counts_ten_train.to(device)
    X_test = counts_ten_test.to(device)

    Z_train = model.encoder(X_train).cpu().numpy()
    Z_test = model.encoder(X_test).cpu().numpy()

print("Training latent space:", Z_train.shape)
print("Test latent space:", Z_test.shape)



latent_cols = [f"latent_{i+1}" for i in range(Z_train.shape[1])]

Z_train_df = pd.DataFrame(
    Z_train,
    index=counts_train_norm.index,
    columns=latent_cols
)

Z_test_df = pd.DataFrame(
    Z_test,
    index=counts_test_norm.index,
    columns=latent_cols
)

Z_train_df.to_csv("./output_data/latent_train.csv")
Z_test_df.to_csv("./output_data/latent_test.csv")
