library(dplyr)

# Extract data

metadata_drug <- read.csv("TCGA-TCIA/Datos/drug_count.csv")

# Dictionary where it establishes to what generic pharmaceutical name will the comercial or the incorrectly written 
# name change to

drug_dictionary <- c(
  "arimidex" = "anastrozole",
  "arimidex (anastrozole)" = "anastrozole",
  "herceptin" = "trastuzumab",
  "trustuzumab" = "trastuzumab",
  "taxol" = "paclitaxel",
  "taxotere" = "docetaxel",
  "cytoxan" = "cyclophosphamide",
  "cyotxan" = "cyclophosphamide",
  "femara" = "letrozole",
  "taxotere/cytoxan" = "docetaxel/cyclophosphamide",
  "aromasin" = "exemestane",
  "avastin" = "bevacizumab"
  
)

# Creates a new column which contains the generic correct names, if the original name doesnt have a coincidence
# with any term from the dictionary it stays the same

metadata_drug2 <- metadata_drug %>%
  mutate(
    clean = recode(pharmaceutical_therapy_drug_name, !!!drug_dictionary)
  )

# Sum of the values in the repeated names on the new column

metadata_drug_limpio <- metadata_drug2 %>%
  group_by(clean) %>%
  summarise(Drug_quantity = sum(Drug_quantity), .groups = "drop") # In the end it eliminates the column created previously
