### Install and load necessary libraries
# If not installed, install required packages
list_of_packages <- c("dplyr", "GGally", "corrplot", "ggplot2", "rgl", "mice", "VIM")
new_packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

# Load libraries
library(dplyr)       # Data manipulation
library(GGally)      # Pair plots
library(corrplot)    # Correlation heatmaps
library(ggplot2)     # Visualization
library(rgl)         # 3D visualization
library(mice)        # Handling missing data
library(VIM)         # Visualizing missing data

### Import the data
KOI_table <- read.csv("KOI_table_2025.03.23_01.50.59.csv", comment.char = "#", header = TRUE, stringsAsFactors = FALSE)

### Data Exploration
head(KOI_table)       # First few rows
names(KOI_table)      # Column names
str(KOI_table)        # Structure of the dataset
dim(KOI_table)        # Dataset dimensions

### Data Cleaning and Preprocessing
# Remove irrelevant columns
columns_to_remove <- c("kepid", "kepoi_name", "kepler_name", "koi_pdisposition")
KOI_table <- KOI_table %>% select(-all_of(columns_to_remove))

# Remove unconfirmed candidate exoplanets
KOI_table <- KOI_table %>% filter(koi_disposition != "CANDIDATE")

# Remove 'koi_score' as it is not relevant
KOI_table <- KOI_table %>% select(-c("koi_score"))

# Encode 'koi_disposition' as binary: 1 (Confirmed), 0 (False Positive)
KOI_table$koi_disposition <- ifelse(KOI_table$koi_disposition == "CONFIRMED", 1, 0)

# Remove non-informative columns
columns_to_remove2 <- c("koi_teq_err1", "koi_teq_err2", "koi_tce_delivname")
KOI_table <- KOI_table %>% select(-all_of(columns_to_remove2))

# Handle missing data
KOI_table_clean <- na.omit(KOI_table)  # Remove rows with NA values

### Data Analysis
# Pairplot for selected features
selected_features <- KOI_table_clean[, c("koi_fpflag_nt", "koi_period", "koi_duration", "koi_teq", "koi_steff", "koi_disposition")]
ggpairs(selected_features, aes(color = as.factor(koi_disposition), alpha = 0.5))

# Correlation matrix
cor_matrix <- cor(KOI_table_clean[, -which(names(KOI_table_clean) == "koi_disposition")], use = "complete.obs")
corrplot(cor_matrix, method = "color", tl.cex = 0.6, tl.srt = 45)

# Boxplots (Standardized Data)
KOI_table_scaled <- KOI_table_clean
KOI_table_scaled[, -which(names(KOI_table_scaled) == "koi_disposition")] <- scale(KOI_table_scaled[, -which(names(KOI_table_scaled) == "koi_disposition")])
boxplot(KOI_table_scaled, las = 2, main = "Standardized Features Boxplot")

# Class distribution
barplot(prop.table(table(KOI_table_clean$koi_disposition)), col = c("red", "blue"), ylim = c(0,1), main = "Class Distribution", ylab = "Proportion", xlab = "KOI Disposition")

### PCA Analysis
pca <- prcomp(KOI_table_clean[, -which(names(KOI_table_clean) == "koi_disposition")], scale. = TRUE)
summary(pca)

# PCA scatter plot with transparency
df_pca <- data.frame(pca$x[, 1:2], label = as.factor(KOI_table_clean$koi_disposition))
ggplot(df_pca, aes(x = PC1, y = PC2, color = label, alpha = 0.5)) +
  geom_point() +
  theme_minimal() +
  ggtitle("PCA: First Two Principal Components")

# PCA for confirmed exoplanets only
KOI_exoplanets <- KOI_table_clean %>% filter(koi_disposition == 1)
df_pca_exo <- data.frame(pca$x[KOI_table_clean$koi_disposition == 1, 1:2])
ggplot(df_pca_exo, aes(x = PC1, y = PC2)) +
  geom_point(color = "blue") +
  xlim(range(df_pca$PC1)) + ylim(range(df_pca$PC2)) +
  theme_minimal() +
  ggtitle("PCA: Confirmed Exoplanets Only")

# 3D PCA plot with transparency
df_pca_3d <- data.frame(pca$x[, 1:3], label = KOI_table_clean$koi_disposition)
plot3d(df_pca_3d$PC1, df_pca_3d$PC2, df_pca_3d$PC3, 
       col = ifelse(df_pca_3d$label == 1, "blue", "red"),
       alpha = 0.5, size = 5, 
       xlab = "PC1", ylab = "PC2", zlab = "PC3", 
       main = "3D PCA Visualization")

# PCA loadings histogram
loadings <- as.data.frame(pca$rotation)  # Loadings matrix
loadings$Feature <- rownames(loadings)   # Feature names
loadings_long <- loadings %>% pivot_longer(cols = starts_with("PC"), names_to = "Principal_Component", values_to = "Loading")

ggplot(loadings_long %>% filter(Principal_Component == "PC1"), aes(x = reorder(Feature, Loading), y = Loading)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +  # Flip axes for readability
  theme_minimal() +
  labs(title = "PCA Loadings for PC1", x = "Features", y = "Loading Value")
