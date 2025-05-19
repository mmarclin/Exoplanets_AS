##########  Study of the outliers ###################################

# Load necessary libraries
library(ggplot2)
library(reshape2)
library(dplyr)       # Data manipulation
library(GGally)      # Pair plots
library(corrplot)    # Correlation heatmaps
library(ggplot2)     # Visualization
library(rgl)         # 3D visualization
library(mice)        # Handling missing data
library(VIM)     

### Import the data ----
KOI_table_full = read.csv("C:/Users/antoi/OneDrive/Bureau/AppliedStats/PROJECT EXOPLANETS/KOI_table_2025.03.23_01.50.59.csv",comment.char = "#",header = TRUE, stringsAsFactors = T)

### Data Cleaning and Preprocessing

# First we need to create the dataset with the target variable and the features 
# the features are all the variables - those who are not useful for the model

columns_to_remove <- c("kepid", "kepoi_name", "kepler_name",
                       "koi_pdisposition","koi_score" ,
                       "koi_teq_err1", "koi_teq_err2", "koi_tce_delivname",
                       "koi_fpflag_nt","koi_fpflag_ss","koi_fpflag_co","koi_fpflag_ec" # variables that indicate how the unit was classified as FP
                       , "ra", "dec" # coordinates of planets in the sky ... not predictive
                       ,"koi_tce_plnt_num") # order in which the object was detected in its star system

KOI_table = KOI_table_full %>% select(-all_of(columns_to_remove))

# Remove unconfirmed candidate exoplanets
KOI_table <- KOI_table %>% filter(KOI_table$koi_disposition != "CANDIDATE")

# Encode 'koi_disposition' as binary: 1 (Confirmed), 0 (False Positive)
KOI_table$koi_disposition <- ifelse(KOI_table$koi_disposition == "CONFIRMED", 1, 0)

KOI_table$koi_disposition = as.factor(KOI_table$koi_disposition)

# Handle missing data : 
# 1st approach : remove rows with NA

KOI_table_clean <- na.omit(KOI_table)  # Remove rows with NA values

summary(KOI_table_clean)


# ----
data = KOI_table_clean
##
numeric_data <- data[sapply(data, is.numeric)]
# Loop through each numeric variable and plot a boxplot
for (feature in names(numeric_data)) {
  p <- ggplot(data, aes_string(y = feature)) +
    geom_boxplot(fill = "lightblue", outlier.colour = "red", outlier.shape = 1) +
    theme_minimal() +
    labs(title = paste("Boxplot of", feature), y = feature)
  
  print(p)
}
# We observe a lot of outliers for almost each numerical feature 
# but real physical measurements 

data[which(data$koi_srad>2),]$koi_disposition
# example srad above a threshold => only FP 

# let's look if the outliers are the same for every features : 
numeric_cols <- names(data)[sapply(data, is.numeric)]
outlier_indices_list <- list()
indexes = 1:7132

for (col in numeric_cols) {
  x <- data[[col]]
  Q1 <- quantile(x, 0.25, na.rm = TRUE)
  Q3 <- quantile(x, 0.75, na.rm = TRUE)
  IQR_val <- Q3 - Q1
  lower_bound <- Q1 - 1.5 * IQR_val
  upper_bound <- Q3 + 1.5 * IQR_val
  outlier_indices_list[[col]] <- which(x < lower_bound | x > upper_bound)
  indexes = setdiff(indexes,which(x < lower_bound | x > upper_bound))
  
}
outlier_indices_list
length(indexes)
# if we remove all outliers we would keep only 2500 datapoints... 
# we shouldn't discard the outliers 


## Outliers in the errors : let's see if they are associated with FP
# if not we should discard them
error_cols <- c("koi_insol_err1", "koi_slogg_err1","koi_steff_err2","koi_steff_err1")  # replace with your actual error columns

for (col in error_cols) {
  p <- ggplot(data, aes_string(x = col)) +
    geom_histogram(bins = 50, fill = "skyblue", color = "black") +
    theme_minimal() +
    labs(title = paste("Distribution of", col), x = col, y = "Count")
  print(p)
}

# check correlation with exoplanet : 
for (col in error_cols) {
  p <- ggplot(data, aes_string(x = "koi_disposition", y = col)) +
    geom_boxplot(fill = "lightgreen") +
    theme_minimal() +
    labs(title = paste(col, "by Label"), x = "koi_disposition", y = col)
  print(p)
}
# we see that higher errors are associated with false positive more often
# we want to find a value to cap the errors, errors extremely high dont bring 
# more info than high errors
sum(data$koi_srad_err1>0.5 & data$koi_disposition==1)
sum(data$koi_srad_err1>0.5 & data$koi_disposition==0)

ggplot(data, aes(x = koi_duration_err1, fill = koi_disposition)) +
  geom_histogram(bins = 100, position = "identity", alpha = 0.5) +
  labs(title = "Distribution of koi_duration_err1 by disposition") +
  theme_minimal()

data_uncapped = data

data$koi_duration_err1 <- pmin(data$koi_duration_err1, 0.5)

# we repeat for each Feature : 

distrib_95 <- list()
features_to_cap = features
# Loop through each feature
for (col in features) {
  val <- quantile(data[[col]], 0.95, na.rm = TRUE)
  
  n_1 <- sum(data[[col]] > val & data$koi_disposition == 1, na.rm = TRUE)
  n_0 <- sum(data[[col]] > val & data$koi_disposition == 0, na.rm = TRUE)
  prop = max(n_1,n_0)/(n_1+n_0)
  if(prop<0.7){
    features_to_cap = setdiff(features_to_cap,col)
  }
  distrib_95[[col]] <- list(
    above_95_class1 = n_1,
    above_95_class0 = n_0,
    prop_class = prop
  )
}

distrib_95
length(features_to_cap)
# for almost every feature, above a certain threshold there is only FP or Only Confirmed
# => We can cap the values of all the features that exhibit this pattern : data_partially_capped
# => we can cap all the features : data_capped
data = data_uncapped

# Computation of the partially capped dataset : 

data_partially_capped = data
for (col in features_to_cap) {
  val <- quantile(data[[col]], 0.95, na.rm = TRUE)
  data_partially_capped[[col]] <- pmin(data[[col]], val)
}


# Computation of the capped dataset : 
data_capped = data
for (col in features) {
  val <- quantile(data[[col]], 0.95, na.rm = TRUE)
  data_capped[[col]] <- pmin(data[[col]], val)
}
#boxplots of the capped features : 

numeric_data <- data[sapply(data, is.numeric)]
# Loop through each numeric variable and plot a boxplot
for (feature in names(numeric_data)) {
  p <- ggplot(data_capped, aes_string(y = feature)) +
    geom_boxplot(fill = "lightblue", outlier.colour = "red", outlier.shape = 1) +
    theme_minimal() +
    labs(title = paste("Boxplot of", feature), y = feature)
  
  print(p)
}


## Other transformations of data to limit influence of outliers : 

# test with log(data) : 

data_log =  log(abs(data[,-c(1)]))

# Pb some have inf:
summary(data_log)

# we can use log(1+x) : 
data_log =  log1p(abs(data[,-c(1)]))

# 
data_log$koi_disposition = data$koi_disposition
numeric_data <- data_log[sapply(data_log, is.numeric)]
# Loop through each numeric variable and plot a boxplot
for (feature in names(numeric_data)) {
  p <- ggplot(data_log, aes_string(y = feature)) +
    geom_boxplot(fill = "lightblue", outlier.colour = "red", outlier.shape = 1) +
    theme_minimal() +
    labs(title = paste("Boxplot of", feature), y = feature)
  
  print(p)
}
# boxplots look better : except koi_slogg => we can cap koi_slogg
f_to_cap = c("koi_slogg","koi_time0bk_err2","koi_time0bk_err1",
             "koi_period_err1","koi_period_err2")
for(col in f_to_cap){
  x = data_log[[col]]
  p5 <- quantile(x, 0.05, na.rm = TRUE)
  p95 <- quantile(x, 0.95, na.rm = TRUE)
  
  data_log[[col]] <- pmin(pmax(x, p5), p95)

}

hist(data_log$koi_period_err1, breaks = 100)
# we that this feature ( period error 1 and 2 ) still have a lot of outliers after 
# log transform and capping, we can simply remove them
# since they are not very discriminative : data_log_f

data_log_f = data_log[,-which(colnames(data_log)==c("koi_period_err2","koi_period_err1"))]


## ----

# We now have the datasets : 

data_log
data_log_f 
data_capped 
data_partially_capped

### Check Correlation of features : 

numeric_data <- data_log[sapply(data_log, is.numeric)]

# Compute correlation matrix
cor_matrix <- cor(numeric_data, use = "pairwise.complete.obs")

# Plot the heatmap
corrplot(cor_matrix, method = "color", tl.cex = 0.6, number.cex = 0.7, type = "upper")
## very high correlation between some features : we need to perform a PCA  

### PCA on new features : 

# first we need to scale our features : 
robust_scale <- function(x) {
  med <- median(x, na.rm = TRUE)
  iqr_val <- IQR(x, na.rm = TRUE)
  return((x - med) / iqr_val)
}
data_scale = data_log[,-c(34)]
data_scale = sapply(data_scale, robust_scale)
#  pca : 
pca_robust <- prcomp(data_scale, scale. = F)
summary(pca_robust)


pca_log <- prcomp(data_log[,-c(34)], scale. = T)
summary(pca_log)


pca_capped <- prcomp(data_capped[,-c(1)], scale. = T)
summary(pca_capped)


# Now we can train some models sensitive to outliers and correlation more safely
