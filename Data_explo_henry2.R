### Install and load necessary libraries
# If not installed, install required packages
list_of_packages <- c("dplyr", "GGally","reshape2" "corrplot", "ggplot2", "rgl", "mice", "VIM")
new_packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

# Load libraries
library(dplyr)       # Data manipulation
library(GGally)      # Pair plots
library(corrplot)    # Correlation heatmaps
library(ggplot2)     # Box plots
library(rgl)         # 3D visualization
library(mice)        # Handling missing data
library(VIM)         # Visualizing missing data

### Import the data
# Ensure the file is in the correct working directory
KOI_table <- read.csv("KOI_table_2025.03.23_01.50.59.csv", comment.char = "#", header = TRUE, stringsAsFactors = FALSE)

### Data Exploration
# Overview of the dataset
head(KOI_table)       # Show first few rows
names(KOI_table)      # Show column names
str(KOI_table)        # Structure of the dataset
dim(KOI_table)        # Dimensions of the dataset
class(KOI_table)      # Type of data object

### Data Cleaning and Preprocessing
# Remove non-relevant features
columns_to_remove <- c("kepid", "kepoi_name", "kepler_name", "koi_pdisposition")
KOI_table <- KOI_table %>% select(-all_of(columns_to_remove))
names(KOI_table)  # Verify removed columns

# Remove unconfirmed candidate exoplanets
KOI_table <- KOI_table %>% filter(koi_disposition != "CANDIDATE")
head(KOI_table)
dim(KOI_table)   # Check new dimensions

# Remove 'koi_score' since we focus on classification
KOI_table <- KOI_table %>% select(-c("koi_score"))
names(KOI_table)

# Encode 'koi_disposition' as binary: 1 (Confirmed), 0 (False Positive)
KOI_table$koi_disposition <- ifelse(KOI_table$koi_disposition == "CONFIRMED", 1, 0)
head(KOI_table)

# Check columns containing strings or NA values
if ("koi_teq_err1" %in% names(KOI_table)) unique(KOI_table$koi_teq_err1)
if ("koi_teq_err2" %in% names(KOI_table)) unique(KOI_table$koi_teq_err2)
if ("koi_tce_delivname" %in% names(KOI_table)) unique(KOI_table$koi_tce_delivname)

# Remove non-informative columns
columns_to_remove2 <- c("koi_teq_err1", "koi_teq_err2", "koi_tce_delivname")
KOI_table <- KOI_table %>% select(-all_of(columns_to_remove2))

# Check for missing values
missing_counts <- colSums(is.na(KOI_table))
print(missing_counts)  # Show columns with missing values

# Visualize missing data
aggr(KOI_table, col=c("navyblue", "red"), numbers=TRUE, sortVars=TRUE, labels=names(KOI_table), cex.axis=.7)

# Handle missing data using multiple imputation (instead of dropping rows)
#KOI_table_imputed <- mice(KOI_table, method = "pmm", m = 5)
#KOI_table_clean <- complete(KOI_table_imputed)
KOI_table_clean <- na.omit(KOI_table_clean)

# Check the final cleaned dataset
dim(KOI_table_clean)
summary(KOI_table_clean)



### Data Analysis

# Overview of the dataset
unique(KOI_table_clean$koi_disposition)  # Check unique values in target variable
names(KOI_table_clean)  # Get feature names
summary(KOI_table_clean)  # Summary statistics
str(KOI_table_clean)  # Structure of the dataset

define_colors <- function(labels) {
  return(rainbow(length(unique(labels))))
}

### Pairplot for selected features


selected_features <- KOI_table_clean[, c("koi_fpflag_nt", "koi_period", "koi_duration", "koi_teq", "koi_steff", "koi_disposition")]

ggpairs(selected_features, aes(color = as.factor(koi_disposition)))  # Corrected color assignment

# Compute correlation matrix
cor_matrix <- cor(KOI_table_clean[, -which(names(KOI_table_clean) == "koi_disposition")], use = "complete.obs")
corrplot(cor_matrix, method = "color", tl.cex = 0.6, tl.srt = 45)

### Boxplots (standardized data for better visualization)
KOI_table_scaled <- KOI_table_clean
KOI_table_scaled[, -which(names(KOI_table_scaled) == "koi_disposition")] <- scale(KOI_table_scaled[, -which(names(KOI_table_scaled) == "koi_disposition")])
boxplot(KOI_table_scaled, las = 2, main = "Standardized Features Boxplot")

### Check for class imbalance
# Create a frequency table
table_koi <- table(KOI_table_clean$koi_disposition)
barplot(prop.table(table_koi),
        col = define_colors(table_koi),
        ylim = c(0,1),
        main = "Class Distribution",
        ylab = "Proportion",
        xlab = "KOI Disposition")

### PCA for Dimensionality Reduction

colSums(is.na(KOI_table_clean))  # Vérifie combien de NA sont présents par colonne
KOI_table_clean <- na.omit(KOI_table_clean)  # Supprime les lignes avec des NA restants

pca <- prcomp(KOI_table_clean[, -which(names(KOI_table_clean) == "koi_disposition")], scale. = TRUE)
summary(pca)


# Select the first 2 principal components
df_pca <- data.frame(pca$x[, 1:2], label = as.factor(KOI_table_clean$koi_disposition))
ggplot(df_pca, aes(x = PC1, y = PC2, color = label)) +
  geom_point() +
  theme_minimal() +
  ggtitle("PCA: First Two Principal Components")

# 3D Scatter Plot
df_pca_3d <- data.frame(pca$x[, 1:3], label = KOI_table_clean$koi_disposition)
plot3d(df_pca_3d$PC1, df_pca_3d$PC2, df_pca_3d$PC3, 
       col = define_colors(df_pca_3d$label)[as.factor(df_pca_3d$label)],
       size = 5, 
       xlab = "PC1", ylab = "PC2", zlab = "PC3", 
       main = "3D PCA Visualization")

# Parallel coordinate plot using first 5 principal components
df_pca_5 <- data.frame(pca$x[, 1:5], label = as.factor(KOI_table_clean$koi_disposition))
ggparcoord(df_pca_5, columns = 1:5, groupColumn = "label", alphaLines = 0.5)




# Calculer la PCA
pca <- prcomp(KOI_table_clean[, -which(names(KOI_table_clean) == "koi_disposition")], scale. = TRUE)

# Extraire les loadings
loadings <- as.data.frame(pca$rotation)  # Matrice des loadings
loadings$Feature <- rownames(loadings)   # Ajouter les noms des variables

# Transformer en format long pour ggplot
loadings_long <- loadings %>%
  pivot_longer(cols = starts_with("PC"), names_to = "Principal_Component", values_to = "Loading")

# Afficher un histogramme des loadings pour PC1
ggplot(loadings_long %>% filter(Principal_Component == "PC1"), aes(x = reorder(Feature, Loading), y = Loading)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +  # Inverser les axes pour une meilleure lisibilité
  theme_minimal() +
  labs(title = "PCA Loadings for PC1",
       x = "Features",
       y = "Loading Value")


