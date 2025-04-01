<<<<<<< HEAD
### KOI_tabl
=======
###
### KOI_table
### Dataset from the NASA : 
### https://exoplanetarchive.ipac.caltech.edu/cgi-bin/TblView/nph-tblView?app=ExoTbls&config=cumulative
###
>>>>>>> 4781138c530b35471f2e8dc19148aa1175cd327e

###
### Import libraries
###

library(gitcreds) # git
library(dplyr) # selection
library(GGally) # pair plots
library(corrplot) # correlation heatmap
library(ggplot2) # box plots
library(rgl) # 3D plot



###
### Import the data 
###

KOI_table = read.csv("KOI_table_2025.03.23_01.50.59.csv",comment.char = "#",header = TRUE, stringsAsFactors = FALSE)

###
### Data exploration
###

head(KOI_table)
names(KOI_table)
str(KOI_table)
dim(KOI_table)
class(KOI_table)

# Remove features we are not interested in
# koi_pdisposition : disposition using only Kepler data
column_rm = c("kepid","kepoi_name","kepler_name","koi_pdisposition")
KOI_table <- KOI_table %>% select(-all_of(column_rm))
names(KOI_table)

# Remove the still unconfirmed candidate exoplanets
# Because we work with only ones whose analyse is finished
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
unique(KOI_table$koi_teq_err2) 
unique(KOI_table$koi_tce_delivname)
column_rm2 = c("koi_teq_err1","koi_teq_err2","koi_tce_delivname")
# We don't have data for the teq err 1 and 2
# tce_delivname is just the delivery type
KOI_table <- KOI_table %>% select(-all_of(column_rm2))

colSums(is.na(KOI_table))
# We can : 
# remove NA
KOI_table_clean = na.omit(KOI_table)
dim(KOI_table_clean)
dim(KOI_table)
# We lose 653 observations if we remove all the NA
# replace NA with specific value (mean, mice)






###
### Data analysing 
###

unique(KOI_table$koi_disposition)
names(KOI_table)
summary(KOI_table)
str(KOI_table)
attach(KOI_table_clean)


### Pairplot for some features
pairs(KOI_table)
selected_features <- KOI_table_clean[, c("koi_fpflag_nt", "koi_period", "koi_duration", "koi_teq", "koi_steff", "koi_disposition")]
ggpairs(selected_features, aes(color = "koi_disposition"))
pairs(selected_features, aes(color="koi_disposition"))

# Compute correlation matrix
cor_matrix <- cor(KOI_table_clean[, -which(names(KOI_table_clean) == "koi_disposition")], use = "complete.obs")
corrplot(cor_matrix, method = "color", tl.cex = 0.6, tl.srt = 45)

### Boxplots
boxplot(KOI_table)
# we should maybe standardize the data

### Check for class imbalance
# Create a frequency table
table_koi <- table(KOI_table_clean[["koi_disposition"]])
table_koi
# Convert to proportions
prop_koi <- prop.table(table_koi)
# Create a bar plot
barplot(prop_koi,
        col = rainbow(length(prop_koi)),  # Assign colors dynamically
        ylim = c(0,1),
        main = "Class Distribution",
        ylab = "Proportion",
        xlab = "KOI Disposition")

### PCA
# We have a lot of features, PCA will help to reduce the dimension and interpret the features
# Compute PCA
pca <- prcomp(KOI_table_clean[, -which(names(KOI_table_clean) == "koi_disposition")], scale. = TRUE)
summary(pca)

# Select the first 2 principal components
df_pca <- data.frame(pca$x[, 1:2], label = koi_disposition)
df_pca
ggplot(df_pca, aes(x = PC1, y = PC2, color = label)) +
  geom_point() +
  theme_minimal()

# 3D Scatter Plot
df_pca <- data.frame(pca$x[, 1:3], label = koi_disposition)
plot3d(df_pca$PC1, df_pca$PC2, df_pca$PC3, 
       col = rainbow(length(unique(df_pca$label)))[as.factor(df_pca$label)],
       size = 5, 
       xlab = "PC1", ylab = "PC2", zlab = "PC3", 
       main = "3D PCA Visualization")

# Select the first 5 principal components
df_pca_5 <- data.frame(pca$x[, 1:5], label = koi_disposition)
# Parallel coordinate plot
ggparcoord(df_pca_5, columns = 1:5, groupColumn = "label", alphaLines = 0.5)

View(KOI_table)

summary(KOI_table)

<<<<<<< HEAD
table(KOI_table$koi_disposition)

detach(KOI_table)
=======


###
### Prediction
###




detach(KOI_table_clean)
>>>>>>> 4781138c530b35471f2e8dc19148aa1175cd327e

library(ggplot2)
library(dplyr)
