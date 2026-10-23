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

counts_train = counts.sample(n=1000)
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

class VAE(nn.Module):
    def __init__(self, input_dim = 5000, latent_dim = 20, hidden_deep = 750, hidden_sup = 3000):
        super().__init__()

        # Encoder

        self.fc0 = nn.Linear(input_dim, hidden_sup) 
        self.fc1 = nn.Linear(hidden_sup, hidden_deep)
        self.fc_mu = nn.Linear(hidden_deep, latent_dim)
        self.fc_logvar = nn.Linear(hidden_deep, latent_dim)

        #Decoder

        self.fc2 = nn.Linear(latent_dim, hidden_deep)
        self.fc3 = nn.Linear(hidden_deep, hidden_sup)
        self.fc4 = nn.Linear(hidden_sup, input_dim)

        self.gelu = nn.GELU()
        self.tanh = nn.Tanh()

    def encode(self, x):
        h0 = self.gelu(self.fc0(x)) 
        h = self.tanh(self.fc1(h0))
        mu = self.fc_mu(h)
        logvar = self.fc_logvar(h)
        return mu, logvar

    def reparameterize(self, mu, logvar):
        std = torch.exp(0.5 * logvar)
        eps = torch.randn_like(std)
        z = mu + eps*std
        return z

    def decoder(self, z):
        h1 = self.tanh(self.fc2(z))
        h2 = self.gelu(self.fc3(h1))
        out = self.fc4(h2)
        return out

    def forward(self, x):
        mu, logvar = self.encode(x)
        z = self.reparameterize(mu, logvar)
        x_reconstructed = self.decoder(z)
        return x_reconstructed, mu, logvar, z

def vae_loss(recon_x, x, mu, logvar, beta=1.0):
    recon_loss = nn.functional.smooth_l1_loss(
        recon_x, x, reduction="none"
    ).sum(dim=1).mean()

    kl_loss = -0.5 * (
        1 + logvar - mu.pow(2) - logvar.exp()
    ).sum(dim=1).mean()

    return recon_loss + beta * kl_loss


model = VAE().to(device)

optimizer = optim.Adam(model.parameters(), lr=1e-3, weight_decay=1e-3)

epochs = 10

def train_batch(batch, model, optimizer):
    model.train()
    optimizer.zero_grad()
    input = batch[0].to(device)
    x_recon, mu, logvar, _ = model(input)
    loss = vae_loss(x_recon, input, mu, logvar)
    loss.backward()
    optimizer.step()
    return loss.item()


@torch.no_grad()
def validate_batch(batch, model):
    model.eval()
    input = batch[0].to(device)
    x_recon, mu, logvar, _ = model(input)
    loss = vae_loss(x_recon, input, mu, logvar)
    return loss.item()


log = Report(epochs)

for epoch in range(epochs):

    N = len(train_loader)
    train_loss = 0
    for batch in train_loader:
        loss = train_batch(batch, model, optimizer)
        train_loss += loss
    train_loss /= N

    val_loss = 0
    for batch in test_loader:
        val_loss += validate_batch(batch, model)
    val_loss /= len(test_loader)

    log.record(pos=epoch + 1, trn_loss=train_loss, val_loss=val_loss, end='\r')
    log.report_avgs(epoch + 1)

log.plot_epochs(log=True)

model.eval()

with torch.no_grad():
    X_train = counts_ten_train.to(device)
    X_test = counts_ten_test.to(device)

    _, mu_train, _, _ = model(X_train)
    _, mu_test, _, _ = model(X_test)
    Z_train = mu_train.cpu().numpy()
    Z_test = mu_test.cpu().numpy()

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

Z_train_df.to_csv("./output_data/latent_train_vae.csv")
Z_test_df.to_csv("./output_data/latent_test_vae.csv")
