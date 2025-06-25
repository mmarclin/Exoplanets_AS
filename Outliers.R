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
library(dplyr)  # make sure dplyr is loaded

### Import the data ----
KOI_table_full = read.csv("C:/Users/antoi/OneDrive/Bureau/AppliedStats/PROJECT EXOPLANETS/KOI_table_2025.03.23_01.50.59.csv",comment.char = "#",header = TRUE, stringsAsFactors = T)

### Data Cleaning and Preprocessing

# First we need to create the dataset with the target variable and the features 
# the features are all the variables - those who are not useful for the model

columns_to_remove <- c("kepid", "kepoi_name", "kepler_name",
                       "koi_pdisposition", "koi_score",
                       "koi_teq_err1", "koi_teq_err2", "koi_tce_delivname",
                       "koi_fpflag_nt", "koi_fpflag_ss", "koi_fpflag_co", "koi_fpflag_ec", # variables that indicate how the unit was classified as FP
                       "ra", "dec", # coordinates of planets in the sky ... not predictive
                       "koi_tce_plnt_num") # order in which the object was detected in its star system

KOI_table <- KOI_table_full[, !(names(KOI_table_full) %in% columns_to_remove)]

# Remove unconfirmed candidate exoplanets
KOI_table <- KOI_table %>% filter(KOI_table$koi_disposition != "CANDIDATE")

# Encode 'koi_disposition' as binary: 1 (Confirmed), 0 (False Positive)
KOI_table$koi_disposition <- ifelse(KOI_table$koi_disposition == "CONFIRMED", 1, 0)

KOI_table$koi_disposition = as.factor(KOI_table$koi_disposition)

# Handle missing data : 
# 1st approach : remove rows with NA

KOI_table_clean <- na.omit(KOI_table)  # Remove rows with NA values

summary(KOI_table_clean)

KOI_table_clean$koi_disposition
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
features =numeric_cols
features_to_cap = features
# Loop through each feature
proptot = NULL
for (col in features) {
  val <- quantile(data[[col]], 0.95, na.rm = TRUE)
  
  n_1 <- sum(data[[col]] > val & data$koi_disposition == 1, na.rm = TRUE)
  n_0 <- sum(data[[col]] > val & data$koi_disposition == 0, na.rm = TRUE)
  prop = max(n_1,n_0)/(n_1+n_0)
  proptot = c(proptot,prop)
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
plot(pca_robust)
summary(pca_robust)


pca_log <- prcomp(data_log[,-c(34)], scale. = T)
summary(pca_log)


pca_capped <- prcomp(data_capped[,-c(1)], scale. = T)
summary(pca_capped)


# Now we can train some models sensitive to outliers and correlation more safely



# Assuming KOI_table_clean is your dataframe

# List of main features and their corresponding error columns
features <- c("koi_period", "koi_time0bk", "koi_impact", "koi_duration", "koi_depth", 
              "koi_prad", "koi_insol", "koi_steff", "koi_slogg", "koi_srad")

upper_errs <- paste0(features, "_err1")
lower_errs <- paste0(features, "_err2")

# Function to compute relative uncertainty per feature for all rows
relative_uncertainty <- function(value, err1, err2) {
  abs_err_mean <- (abs(err1) + abs(err2)) / 2
  rel_unc <- abs_err_mean / abs(value)
  return(rel_unc)
}

# Initialize a matrix to store relative uncertainties for each feature
rel_unc_mat <- matrix(NA, nrow = nrow(KOI_table_clean), ncol = length(features))
colnames(rel_unc_mat) <- features

# Compute relative uncertainties for each feature
for (i in seq_along(features)) {
  val <- KOI_table_clean[[features[i]]]
  err1 <- KOI_table_clean[[upper_errs[i]]]
  err2 <- KOI_table_clean[[lower_errs[i]]]
  
  # Avoid division by zero or NA by replacing zeros with small value or NA handling
  val[val == 0] <- NA
  
  rel_unc_mat[, i] <- relative_uncertainty(val, err1, err2)
}

# Aggregate uncertainty across features per observation (mean of available relative uncertainties)
total_uncertainty <- apply(rel_unc_mat, 1, function(x) mean(x, na.rm = TRUE))

# Compute weights: higher weight for lower uncertainty
weights <- 1 / (1 + total_uncertainty)
# Normalize weights between 0 and 1
weights_norm <- weights / max(weights, na.rm = TRUE)

#### ESI : 

data_new = KOI_table_clean
# Reference Earth values
earth_radius <- 1      # in Earth radii
earth_temp <- 288      # in Kelvin
earth_insol <- 1  # Earth receives 1 Earth flux

# Similarity function
param_similarity <- function(x, x_earth) {
  1 - abs((x - x_earth) / (x + x_earth))
}


# Adding insolation to ESI
w_insol <- 1/3
w_radius <- 1/3
w_temp <- 1/3

data_new$ESI <- with(KOI_table_clean,
                            ifelse(!is.na(koi_prad) & !is.na(koi_teq) & !is.na(koi_insol),
                                   (param_similarity(koi_prad, earth_radius)^w_radius) *
                                     (param_similarity(koi_teq, earth_temp)^w_temp) *
                                     (param_similarity(koi_insol, earth_insol)^w_insol),
                                   NA))

summary(data_new$ESI)

KOI_table_clean$ESI = data_new$ESI
### Plot of the ESI : 
KOI_high_esi <- KOI_table_clean %>% 
  filter(ESI > 0.5)

# Plot
ggplot(KOI_high_esi, aes(x = koi_prad, y = koi_teq, color = ESI)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_gradient(low = "yellow", high = "darkblue") +
  labs(title = "Planets with ESI > 0.5",
       x = "Planet Radius (R_Earth)",
       y = "Equilibrium Temperature (K)",
       color = "ESI") +
  theme_minimal()

# Plot density of ESI by disposition
ggplot(KOI_table_clean, aes(x = ESI, fill = koi_disposition)) +
  geom_density(alpha = 0.5) +
  labs(title = "Density of ESI by KOI Disposition",
       x = "Earth Similarity Index (ESI)",
       y = "Density",
       fill = "Disposition") +
  theme_minimal()




# Load the dataset
TESS_table <- read.csv("C:/Users/antoi/OneDrive/Bureau/AppliedStats/PROJECT EXOPLANETS/TOI_2025.06.24_02.38.34.csv", 
                       comment.char = "#", header = TRUE, stringsAsFactors = FALSE)

# Rename TESS columns to match KOI feature names
TESS_table <- TESS_table %>%
  rename(
    koi_disposition = tfopwg_disp,
    koi_period = pl_orbper,
    koi_period_err1 = pl_orbpererr1,
    koi_period_err2 = pl_orbpererr2,
    koi_time0bk = pl_tranmid,
    koi_time0bk_err1 = pl_tranmiderr1,
    koi_time0bk_err2 = pl_tranmiderr2,
    koi_duration = pl_trandurh,
    koi_duration_err1 = pl_trandurherr1,
    koi_duration_err2 = pl_trandurherr2,
    koi_depth = pl_trandep,
    koi_depth_err1 = pl_trandeperr1,
    koi_depth_err2 = pl_trandeperr2,
    koi_prad = pl_rade,
    koi_prad_err1 = pl_radeerr1,
    koi_prad_err2 = pl_radeerr2,
    koi_teq = pl_eqt,
    koi_insol = pl_insol,
    koi_insol_err1 = pl_insolerr1,
    koi_insol_err2 = pl_insolerr2,
    koi_steff = st_teff,
    koi_steff_err1 = st_tefferr1,
    koi_steff_err2 = st_tefferr2,
    koi_slogg = st_logg,
    koi_slogg_err1 = st_loggerr1,
    koi_slogg_err2 = st_loggerr2,
    koi_srad = st_rad,
    koi_srad_err1 = st_raderr1,
    koi_srad_err2 = st_raderr2,
    koi_kepmag = st_tmag
  )

# Add missing variables (fill with zero)
TESS_table$koi_impact <- 0
TESS_table$koi_impact_err1 <- 0
TESS_table$koi_impact_err2 <- 0
TESS_table$koi_insol_err2  <- 0
TESS_table$koi_insol_err1  <- 0

# Weights (you can adjust if needed)
w_insol <- 1/3
w_radius <- 1/3
w_temp <- 1/3

TESS_table$ESI <- with(TESS_table,
                             ifelse(!is.na(koi_prad) & !is.na(koi_teq) & !is.na(koi_insol),
                                    (param_similarity(koi_prad, earth_radius)^w_radius) *
                                      (param_similarity(koi_teq, earth_temp)^w_temp) *
                                      (param_similarity(koi_insol, earth_insol)^w_insol),
                                    NA)
)


# Select only the features you use in the model
model_features <- c(
  "koi_disposition", "koi_period", "koi_period_err1", "koi_period_err2", 
  "koi_time0bk", "koi_time0bk_err1", "koi_time0bk_err2", 
  "koi_impact", "koi_impact_err1", "koi_impact_err2", 
  "koi_duration", "koi_duration_err1", "koi_duration_err2", 
  "koi_depth", "koi_depth_err1", "koi_depth_err2", 
  "koi_prad", "koi_prad_err1", "koi_prad_err2", 
  "koi_teq", "koi_insol", "koi_insol_err1", "koi_insol_err2", 
  "koi_model_snr", "koi_steff", "koi_steff_err1", "koi_steff_err2", 
  "koi_slogg", "koi_slogg_err1", "koi_slogg_err2", 
  "koi_srad", "koi_srad_err1", "koi_srad_err2", 
  "koi_kepmag", "ESI"
)

TESS_table_model <- TESS_table  %>% dplyr::select(all_of(model_features))
summary(TESS_table_model)


##### Change the output
# Map TESS disposition to KOI binary labels (1 = confirmed, 0 = false positive)
TESS_table_model <- TESS_table_model %>%
  mutate(koi_disposition_binary = case_when(
    koi_disposition %in% c("CP", "PC") ~ 1,
    koi_disposition %in% c("FP") ~ 0,
    TRUE ~ NA_real_ # For APC, FA, KP etc. assign NA (to filter out later)
  ))

# Optionally: remove rows where the target is NA
TESS_table_model <- TESS_table_model %>% filter(!is.na(koi_disposition_binary))

# Drop the original 'koi_disposition' and replace with the new binary version
TESS_table_model <- TESS_table_model %>%
  dplyr::select(-koi_disposition) %>%
  rename(koi_disposition = koi_disposition_binary)
#####

TESS_clean = na.omit(TESS_table_model)
TESS_log = log1p(abs(TESS_clean[,-c(35)]))




# Let's start with models using the dataset after log transformation of features : 
data_model = data_log

## Train-Test split : 
set.seed(42)
train_index <- createDataPartition(data_model$koi_disposition, p = 0.8, list = FALSE)
train_data <- data_model[train_index, ]
test_data  <- data_model[-train_index, ]
## Logistic Regression : 
log_model <- glm(koi_disposition ~ ., data = train_data, family = binomial)
summary(log_model)
# Predict on test set (probabilities then classes)
log_probs <- predict(log_model, newdata = test_data, type = "response")
log_pred <- ifelse(log_probs > 0.5, 1, 0)
# Confusion matrix
table(Predicted = log_pred, Actual = test_data$koi_disposition)


## Predict on the tess dataset : 
log_probs_tess <- predict(log_model, newdata = TESS_log, type = "response")
log_pred_tess <- ifelse(log_probs_tess > 0.5, 1, 0)
table(Predicted = log_pred_tess, Actual = TESS_clean$koi_disposition)




### Clustering based on the built up dissimilarity matrix : 
# Assuming KOI_table_clean contains the columns: koi_prad, koi_teq, koi_insol

# Similarity function (between x and y, where y replaces Earth in your original formula)


#### PCA on KOI_table.log_f
summary(KOI_table.log_f)

boxplot(KOI_table.log_f, las = 2, col = 'gold')

# Note: PCA is not about the mean, it is about the variability
boxplot(scale(x = KOI_table.log_f[, -which(names(KOI_table.log_f) == "koi_disposition")], center = T, scale = T), las = 2, col = 'gold')
KOI_scaled = scale(x = KOI_table.log_f[, -which(names(KOI_table.log_f) == "koi_disposition")], center = T, scale = T)
pc.koi <- princomp(KOI_scaled, scores = T)
summary(pc.koi)
# Loadings
load.tour <- pc.koi$loadings
load.tour
load.tour[, 1:8]

#" 
eigenvalues <- pc.koi$sdev^2
barplot(eigenvalues, 
        main = "Scree Plot", 
        xlab = "Principal Components", 
        ylab = "Eigenvalue", 
        col = "skyblue", 
        names.arg = paste0("PC", 1:length(eigenvalues)))
#
cum_var <- cumsum(pc.koi$sdev^2 / sum(pc.koi$sdev^2))
plot(cum_var, type = 'b', pch = 19, col = 'darkgreen',
     xlab = "Number of Principal Components", 
     ylab = "Cumulative Proportion of Variance Explained",
     main = "Cumulative Variance Explained")
abline(h = 0.95, col = "red", lty = 2)  # 95% => we keep 13 components 

# Extract PCA-transformed dataset
KOI_pca_scores <- pc.koi$scores  # Data in the new PCA space
# Optionally bind the disposition column back
KOI_pca <- data.frame(KOI_pca_scores[,1:13], koi_disposition = KOI_table.log_f$koi_disposition)
# View the PCA-transformed dataset
head(KOI_pca)

