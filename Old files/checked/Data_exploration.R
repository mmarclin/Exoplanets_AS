###
### KOI_table
### Dataset from the NASA : 
### https://exoplanetarchive.ipac.caltech.edu/cgi-bin/TblView/nph-tblView?app=ExoTbls&config=cumulative
###

###
### Import libraries
###

library(gitcreds) # git
library(dplyr) # selection
library(GGally) # pair plots
library(corrplot) # correlation heatmap
library(ggplot2) # box plots
library(rgl) # 3D plot

library(caret) # split data
library(class) # KNN

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

unique(KOI_table_clean$koi_disposition)
names(KOI_table_clean)
summary(KOI_table_clean)
str(KOI_table)
attach(KOI_table_clean)
list_features = setdiff(colnames(KOI_table_clean),c("koi_disposition"))

### Outliers

outlier_indices_list <- list()  # Initialize empty list to store indices

for (feature_name in names(KOI_table_clean)) {
  feature <- KOI_table_clean[[feature_name]]  # Extract column
  
  if (is.numeric(feature)) {
    # Get boxplot stats
    stats <- boxplot.stats(feature)
    
    # Identify indices of outliers
    outlier_indices <- which(feature %in% stats$out)
    
    # Store in list with feature name as key
    outlier_indices_list[[feature_name]] <- outlier_indices
    
    # Plot boxplot with labels
    boxplot(feature,
            main = paste("Boxplot of", feature_name),
            xlab = feature_name)
  }
}

# Filter out empty entries (features with no outliers)
non_empty_lists <- Filter(length, outlier_indices_list)
# Get the union of all outlier indices
all_outliers <- unique(unlist(non_empty_lists))
all_outliers
dim(all_outliers)
# We lose too much information by removing the outliers
# We can do transformations 

# We can't use log or box-cox transformation because we also have negative values
# However, we can use Yeo-Johnson transformation
library(e1071)
not_skewness_features <- list()
for (feature in KOI_table_clean) {
  hist(feature)
}

for 
skewness(df$feature1)


# Per class
for (feature in KOI_table_clean){
  boxplot(feature ~ koi_disposition, data = KOI_table_clean, main = "Boxplot by Class")
}
dim(koi_disposition)


### Pairplot for some features
pairs(KOI_table)
selected_features <- KOI_table_clean[, c("koi_fpflag_nt", "koi_period", "koi_duration", "koi_teq", "koi_steff", "koi_disposition")]
ggpairs(selected_features, aes(color = "koi_disposition"))
pairs(selected_features, aes(color="koi_disposition"))


### Density 
for (var in list_features) {
  print(
    ggplot(KOI_table_clean, aes(x = !!sym(var), fill = factor(koi_disposition))) +
      geom_density(alpha = 0.5) +
      scale_fill_manual(values = c("steelblue4", "indianred", "darkgreen")) +
      theme_minimal() +
      ggtitle(paste("Density Plot of", var))
  )
}


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

# Cumulative variance
# Variance explained
var_explained <- pca$sdev^2 / sum(pca$sdev^2)
cumvar <- cumsum(var_explained)
plot(cumvar, 
     xlab = "Number of Principal Components", 
     ylab = "Cumulative Proportion of Variance Explained",
     type = "b", pch = 19, col = "blue",
     main = "Cumulative Variance Explained by PCs")
abline(h = 0.9, col = "red", lty = 2)  # optional 90% guide


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





###
### Prediction
###

### Cross validation
set.seed(123)  # For reproducibility
train_index <- createDataPartition(koi_disposition, p = 0.8, list = FALSE)  # 80% training

train_set <- KOI_table_clean[train_index, ]
val_set <- KOI_table_clean[-train_index, ]

head(train_set)
dim(train_set)
dim(val_set)

### LDA

KOI_table_clean.knn <- knn(train = , cl = koi_disposition, k = 3, prob = T)


### Logistic regression

glm.fit <- glm(koi_disposition~., data = train_set, family=binomial)
summary(glm.fit)

# Prediction on validation set
val_probs <- predict(glm.fit, newdata = val_set, type = "response")
# Class prediction
val_pred <- ifelse(val_probs > 0.5, 1, 0)
# Confusion matrix
table(Predicted = val_pred, Actual = val_set$koi_disposition)
# Accuracy
mean(val_pred == val_set$koi_disposition)


### KNN




### SVM


detach(KOI_table_clean)


