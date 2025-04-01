### KOI_tabl

library(gitcreds)
library(dplyr)

# Data importation
KOI_table = read.csv("KOI_table_2025.03.23_01.50.59.csv",comment.char = "#",header = TRUE, stringsAsFactors = FALSE)

# Data exploration
head(KOI_table)
names(KOI_table)
str(KOI_table)
dim(KOI_table)
class(KOI_table)

# Remove features we are not interested in
column_rm = c("kepid","kepoi_name","kepler_name","koi_pdisposition")
KOI_table <- KOI_table %>% select(-all_of(column_rm))
names(KOI_table)

# Remove the still unconfirmed candidate exoplanets
KOI_table <- KOI_table %>% filter(koi_disposition != "CANDIDATE")
head(KOI_table)
dim(KOI_table)

# We can work on the score or on the class of the exoplanet (confirmed and false positive)
# Let's choose the class
KOI_table <- KOI_table %>% select(-c("koi_score"))
names(KOI_table)
# Binary encoding : 0 if false positive, 1 if confirmed exoplanet
KOI_table$koi_disposition <- ifelse(KOI_table$koi_disposition == "CONFIRMED", 1, 0)
head(KOI_table)

# We have some columns with NA and some with strings we have to check
# COLUMN koi_depth_err1: Transit Depth Upper Unc. [ppm]
# COLUMN koi_depth_err2: Transit Depth Lower Unc. [ppm]
# COLUMN koi_tce_delivname: TCE Delivery
unique(KOI_table$koi_teq_err1)
unique(KOI_table$koi_teq_err2) #may be useful : to check
unique(KOI_table$koi_tce_delivname)
column_rm2 = c("koi_teq_err1","kepoi_name","kepler_name","koi_pdisposition")

# Data exploration
unique(KOI_table$koi_disposition)
summary(KOI_table)
pairs(KOI_table)
attach(KOI_table)

View(KOI_table)

summary(KOI_table)

table(KOI_table$koi_disposition)

detach(KOI_table)

library(ggplot2)
library(dplyr)

# Filtrer les données avec période < 100
subset_data <- KOI_table %>% filter(koi_period < 100)

# Créer un histogramme avec deux couleurs
ggplot(subset_data, aes(x = koi_period, fill = factor(koi_disposition))) +
  geom_histogram(binwidth = 5, position = "identity", alpha = 0.6, color = "black") +
  scale_fill_manual(values = c("red", "blue"), labels = c("Non confirmée", "Confirmée")) +
  labs(title = "Distribution des périodes (< 100 jours)", x = "Période (jours)", y = "Nombre d'observations", fill = "Statut") +
  theme_minimal()

# Sélectionner les 6 premières variables numériques (exclure koi_disposition)
numeric_features <- subset_data %>% select(where(is.numeric)) %>% select(-koi_disposition) %>% select(1:6)

# Boucle pour générer les histogrammes pour ces variables
for (feature in colnames(numeric_features)) {
  p <- ggplot(subset_data, aes(x = .data[[feature]], fill = factor(koi_disposition))) +
    geom_histogram(binwidth = 5, position = "identity", alpha = 0.6, color = "black") +
    scale_fill_manual(values = c("red", "blue"), labels = c("Non confirmée", "Confirmée")) +
    labs(title = paste("Distribution de", feature), x = feature, y = "Nombre d'observations", fill = "Statut") +
    theme_minimal()
  
  print(p)  # Affiche chaque histogramme
}

