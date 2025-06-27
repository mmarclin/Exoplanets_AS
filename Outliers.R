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


### Habitability prediction using HWC dataset : 

hwc = read.csv("C:/Users/antoi/OneDrive/Bureau/AppliedStats/PROJECT EXOPLANETS/hwc.csv")
head(hwc)
# Define the features to keep
selected_features <- c(
  # Planet
  "P_RADIUS", "P_PERIOD", "P_TEMP_EQUIL", "P_FLUX",
  
  # Star
  "S_TEMPERATURE", "S_RADIUS", "S_LOG_G",
  
  # Target
  "P_HABITABLE"
)

# Keep only those columns in the dataset
hwc <- hwc[, selected_features]
hwc$P_HABITABLE = ifelse(hwc$P_HABITABLE ==0 , 0 ,1 )
hwc$P_HABITABLE = as.factor(hwc$P_HABITABLE)
levels(hwc$P_HABITABLE )
table(hwc$P_HABITABLE)
# change the names : 
colnames(hwc) <- c(
  "koi_prad", "koi_period", "koi_teq", "koi_insol",
  "koi_steff", "koi_srad", "koi_slogg",
  "koi_disposition"
)
hwc = na.omit(hwc)
nrow(hwc)
hwc = hwc[-4161,]
# check the scale :

ncol(hwc)
boxplot(scale(x = hwc[which(hwc$koi_disposition!=0),-9], center = T, scale = T), las = 2, col = 'gold')

boxplot(scale(x = hwc[,-9], center = T, scale = T), las = 2, col = 'gold')

boxplot(hwc[which(hwc$koi_disposition!=0),], las = 2, col = 'gold',scale=T,center = T)
hwc$koi_period

# List of features to plot
features <- c("koi_prad", "koi_period", "koi_teq", "koi_insol", 
              "koi_steff", "koi_srad", "koi_slogg")
KOI_table_exo = KOI_table[which(KOI_table$koi_disposition==1),]
hwc_exo = hwc[which(hwc$koi_disposition==1),]

# Loop through each feature
for (feature in features) {
  # Set up the plotting window: 1 row, 2 columns
  par(mfrow = c(1, 2))
  
  # Boxplot for HWC
  boxplot(hwc_exo[[feature]],
          main = paste("HWC -", feature),
          ylab = feature,
          col = "lightblue",
          outline = TRUE)
  
  # Boxplot for KOI
  boxplot(KOI_table_exo[[feature]],
          main = paste("KOI -", feature),
          ylab = feature,
          col = "lightgreen",
          outline = TRUE)
  cat(mean(KOI_table_exo[[feature]]))
  cat("--")
  cat(mean(hwc_exo[[feature]]))
  
  # Pause between plots so the user can see each one
  readline(prompt = "Press [Enter] to show next feature...")
}

set.seed(123)

# Assume hwc_data is your dataset and
# target is 'P_HABITABLE' (0/1)
hwc_data = hwc
# Split data into train and test (e.g., 70-30 split)
train_index <- createDataPartition(hwc_data$koi_disposition, p = 0.7, list = FALSE)
train_data <- hwc_data[train_index, ]
test_data <- hwc_data[-train_index, ]

# Undersample majority class in training set
minority_class <- train_data %>% filter(koi_disposition == 1)
majority_class <- train_data %>% filter(koi_disposition == 0)

set.seed(123)
majority_down <- majority_class %>% sample_n(nrow(minority_class))

train_balanced <- bind_rows(minority_class, majority_down)

# Separate features and labels for XGBoost
train_matrix <- train_balanced %>% select(-koi_disposition) %>% as.matrix()
train_label <- train_balanced$koi_disposition

test_matrix <- test_data %>% select(-koi_disposition) %>% as.matrix()
test_label <- test_data$koi_disposition

# 1) Random Forest

rf_model <- randomForest(as.factor(koi_disposition) ~ ., data = train_balanced)
rf_preds <- predict(rf_model, newdata = test_data)

confusionMatrix(rf_preds, as.factor(test_label))

# 2) XGBoost

dtrain <- xgb.DMatrix(data = train_matrix, label = as.numeric(train_label)-1)
dtest <- xgb.DMatrix(data = test_matrix, label = as.numeric(test_label)-1)

params <- list(
  objective = "binary:logistic",
  eval_metric = "logloss",
  max_depth = 6,
  eta = 0.1
)

xgb_model <- xgb.train(params, dtrain, nrounds = 100)

xgb_preds_prob <- predict(xgb_model, dtest)
xgb_preds <- ifelse(xgb_preds_prob > 0.5, 1, 0)

confusionMatrix(factor(xgb_preds), factor(test_label))

# XGBoost performs well ... 
# let's use it on the KOI dataset : 
# Select features used in training
features <- c("koi_prad", "koi_period", "koi_teq", "koi_insol", 
              "koi_steff", "koi_srad", "koi_slogg")
setdiff(features, colnames(KOI_table_exo))

# Prepare KOI data
#KOI_table_exo = KOI_table
KOI_pred_matrix <- KOI_table_exo[, features] %>% as.matrix()

# Create DMatrix for XGBoost
KOI_dmatrix <- xgb.DMatrix(data = KOI_pred_matrix)
colnames(KOI_dmatrix)
colnames(dtrain)
# Predict using trained model
KOI_preds_prob <- predict(xgb_model, KOI_dmatrix)
KOI_preds <- ifelse(KOI_preds_prob > 0.5, 1, 0)

# Attach predictions to the KOI data
KOI_table_exo$Predicted_Habitable <- KOI_preds

# Optional: View predicted habitables
table(KOI_table_exo$Predicted_Habitable)


ggplot(KOI_table_exo, aes(x = koi_prad, y = koi_teq, color = factor(Predicted_Habitable))) +
  geom_point(alpha = 0.7) +
  labs(title = "Predicted Habitability on KOI Dataset", color = "Predicted\nHabitable") +
  theme_minimal()


ggplot(KOI_table_exo, aes(x = koi_teq, y = ESI, color = factor(Predicted_Habitable))) +
  geom_point(alpha = 0.7) +
  labs(title = "Predicted Habitability on KOI Dataset", color = "Predicted\nHabitable") +
  theme_minimal()


### Nice plots : 
# Assuming you have a vector KOI_preds with values 0 (not habitable) or 1 (habitable)
library(ggplot2)

KOI_pred_df <- data.frame(Habitability = factor(KOI_preds, labels = c("Not Habitable", "Habitable")))

ggplot(KOI_pred_df, aes(x = Habitability, fill = Habitability)) +
  geom_bar() +
  scale_fill_manual(values = c("lightcoral", "lightblue")) +
  theme_minimal() +
  labs(title = "Predicted Habitability (KOI Dataset)", x = "", y = "Count") +
  theme(legend.position = "none")

# Features to plot
features <- c("koi_prad", "koi_period", "koi_teq", "koi_insol", 
              "koi_steff", "koi_srad", "koi_slogg")

# Add predicted labels to KOI table if not already done
KOI_table_exo$habitability_pred <- factor(KOI_preds, labels = c("Not Habitable", "Habitable"))

# Plot boxplots
library(ggplot2)
for (feature in features) {
  p <- ggplot(KOI_table_exo, aes_string(x = "habitability_pred", y = feature, fill = "habitability_pred")) +
    geom_boxplot(outlier.colour = "red", alpha = 0.7) +
    scale_fill_manual(values = c("lightcoral", "lightblue")) +
    theme_minimal() +
    labs(title = paste("Boxplot of", feature, "by Habitability"), x = "", y = feature) +
    theme(legend.position = "none")
  
  print(p)
}
#### 
ggplot(KOI_table, aes(x = factor(koi_disposition))) +
  geom_bar(fill = c("skyblue","red")) +
  labs(title = "Class Distribution (Confirmed vs False Positives)", x = "Disposition", y = "Count") +
  theme_minimal()



#### Density plots in dataExploration : # Plot density of ESI by disposition
ggplot(KOI_table, aes(x = ESI, fill = as.factor(koi_disposition))) +
  geom_density(alpha = 0.5) +
  labs(title = "Density of ESI by KOI Disposition",
       x = "Earth Similarity Index (ESI)",
       y = "Density",
       fill = "Disposition") +
  theme_minimal()



# Get feature names (excluding the target/class variable)
features <- setdiff(names(KOI_table), "koi_disposition")

# Loop through features and plot densities
for (feature in features) {
  if (is.numeric(KOI_table[[feature]])) {
    p <- ggplot(KOI_table, aes_string(x = feature, fill = "as.factor(koi_disposition)")) +
      geom_density(alpha = 0.5) +
      labs(
        title = paste("Density of", feature, "by KOI Disposition"),
        x = feature,
        y = "Density",
        fill = "Disposition"
      ) +
      theme_minimal()
    
    print(p)  # display the plot
  }
}

library(ggplot2)
library(patchwork)  # install.packages("patchwork") if needed

p1 <- ggplot(KOI_table, aes(x = ESI, fill = as.factor(koi_disposition))) +
  geom_density(alpha = 0.5) +
  labs(title = "Density of ESI by KOI Disposition",
       x = "Earth Similarity Index (ESI)",
       y = "Density",
       fill = "Disposition") +
  theme_minimal()

p2 <- ggplot(KOI_table, aes(x = koi_slogg, fill = as.factor(koi_disposition))) +
  geom_density(alpha = 0.5) +
  labs(title = "Density of Surface Gravity by KOI Disposition",
       x = "Log Surface Gravity (koi_slogg)",
       y = "Density",
       fill = "Disposition") +
  theme_minimal()

# Combine side-by-side
p1 + p2
#### 3 dplot : 

# Load libraries
library(MASS)
library(plotly)

# Clean and filter dataset for relevant features and remove NAs
data_3d <- na.omit(KOI_table[, c("koi_slogg", "koi_prad", "koi_disposition")])

# Ensure class label is a factor
data_3d$koi_disposition <- as.factor(data_3d$koi_disposition)

# Create kernel density estimations for each class
dens_list <- lapply(levels(data_3d$koi_disposition), function(class_label) {
  class_data <- subset(data_3d, koi_disposition == class_label)
  kde2d(class_data$koi_slogg, class_data$koi_prad, n = 80)
})

# Assign colors for each class
colors <- c("blue", "red")
names(colors) <- levels(data_3d$koi_disposition)

# Create individual surfaces
plot <- plot_ly()
for (i in seq_along(dens_list)) {
  plot <- add_surface(plot,
                      x = dens_list[[i]]$x,
                      y = dens_list[[i]]$y,
                      z = dens_list[[i]]$z,
                      showscale = FALSE,
                      opacity = 0.6,
                      name = paste("Class", levels(data_3d$koi_disposition)[i]),
                      surfacecolor = matrix(rep(i, length(dens_list[[i]]$x) * length(dens_list[[i]]$y)), 
                                            nrow = length(dens_list[[i]]$x)),
                      colorscale = list(c(0, 1), c(colors[i], colors[i])))
}

# Layout
plot <- plot %>%
  layout(title = "3D Density of koi_slogg vs ESI by Class",
         scene = list(
           xaxis = list(title = "koi_slogg"),
           yaxis = list(title = "ESI"),
           zaxis = list(title = "Density")
         ))

# Show plot
plot


####### PLot stellar MAP : 
KOI_table_full = read.csv("C:/Users/antoi/OneDrive/Bureau/AppliedStats/PROJECT EXOPLANETS/KOI_table_2025.03.23_01.50.59.csv",comment.char = "#",header = TRUE, stringsAsFactors = T)
sum(KOI_table_full$kepler_name == "Kepler-37 d")
### Data Cleaning and Preprocessing

# First we need to create the dataset with the target variable and the features 
# the features are all the variables - those who are not useful for the model

columns_to_remove <- c("kepid", "kepler_name",
                       "koi_pdisposition", "koi_score",
                       "koi_teq_err1", "koi_teq_err2", "koi_tce_delivname",
                       "koi_fpflag_nt", "koi_fpflag_ss", "koi_fpflag_co", "koi_fpflag_ec", # variables that indicate how the unit was classified as FP
                       "koi_tce_plnt_num") # order in which the object was detected in its star system

KOI_table <- KOI_table_full[, !(names(KOI_table_full) %in% columns_to_remove)]

KOI_table <- KOI_table %>% filter(KOI_table$koi_disposition != "CANDIDATE")

library(ggplot2)

# Replace 'koi_disposition' with the name of your class variable if different
# Make sure koi_disposition is a factor with meaningful labels
KOI_table$koi_disposition <- factor(KOI_table$koi_disposition,
                                    levels = c("FALSE POSITIVE", "CONFIRMED"),
                                    labels = c("False Positive", "Confirmed"))

# Plot RA vs Dec
ggplot(KOI_table, aes(x = ra, y = dec, color = koi_disposition)) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_color_manual(values = c("blue",  "red")) +
  labs(title = "Sky Map of KOIs by Disposition",
       x = "Right Ascension (degrees)",
       y = "Declination (degrees)",
       color = "Disposition") +
  theme_minimal()

########## New dataset for clustering : 
KOI_names <- toupper(KOI_table_full$kepler_name)

ps = read.csv("C:/Users/antoi/OneDrive/Bureau/AppliedStats/PROJECT EXOPLANETS/PS_2025.06.26_09.55.40.csv",comment.char = "#",header = TRUE, stringsAsFactors = FALSE)
kepnames = na.omit(KOI_table_full$kepler_name)
head(ps)
dim(ps)
colnames(ps)
ps$pl_name_upper <- toupper(ps$pl_name)

##### Step 2: Filter the confirmed dataset
#exo_koi_overlap <- ps %>%
#  filter(pl_name_upper %in% KOI_names)

# Step 3: Optionally drop the helper column
exo_koi_overlap = ps
exo_koi_overlap <- exo_koi_overlap %>% select(-pl_name_upper)
exo_koi_overlap = exo_koi_overlap[which(exo_koi_overlap$default_flag==1),]
# View the result
nrow(exo_koi_overlap)
head(exo_koi_overlap)
dim(exo_koi_overlap)
#### Only the KOI exoplanets kept  : 
## what features to use : 
# Required libraries
library(tidyverse)
library(GGally)

# Select features
selected_vars <- exo_koi_overlap %>%
  select(pl_rade, pl_bmasse, pl_orbper, pl_orbeccen,
         pl_eqt, pl_insol, st_teff, st_rad, st_mass) 
dim(selected_vars)
# Pairplot
ggpairs(selected_vars)

# Hierarchical clustering : 

# Kmeans : 
# Load libraries
library(factoextra)  # for visualization
sum(is.na(selected_vars))
ps_na = na.omit(selected_vars)
dim(ps_na)
# 1. Scale the variables (important for KMeans)
scaled_vars <- scale(ps_na)

# 2. Determine optimal number of clusters (optional but useful)
fviz_nbclust(scaled_vars, kmeans, method = "wss") + 
  labs(title = "Elbow Method for Optimal k")

# 3. Run KMeans (e.g. with 3 clusters)
set.seed(123)  # for reproducibility
kmeans_result <- kmeans(scaled_vars, centers = 4, nstart = 25)

# 4. Add cluster labels to your dataset
ps_na$cluster <- as.factor(kmeans_result$cluster)

# 5. Visualize clusters (using PCA for 2D)
fviz_cluster(kmeans_result, data = scaled_vars, 
             geom = "point", ellipse.type = "norm", 
             main = "KMeans Clustering")

ggplot(ps_na, aes(x = pl_rade, y = pl_bmasse, color = cluster)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_brewer(palette = "Set2") +
  labs(title = "KMeans Clusters by Radius and Mass",
       x = "Planet Radius (Earth radii)",
       y = "Planet Mass (Earth masses)",
       color = "Cluster") +
  theme_minimal()


scale_x_log10() + scale_y_log10()
centers <- as.data.frame(kmeans_result$centers)
centers$pl_rade <- centers$pl_rade * attr(scaled_vars, "scaled:scale")[1] + attr(scaled_vars, "scaled:center")[1]
centers$pl_bmasse <- centers$pl_bmasse * attr(scaled_vars, "scaled:scale")[2] + attr(scaled_vars, "scaled:center")[2]

ggplot(ps_na, aes(x = pl_rade, y = pl_bmasse, color = cluster)) +
  geom_point(alpha = 0.7) +
  geom_point(data = centers, aes(x = pl_rade, y = pl_bmasse), color = "black", size = 4, shape = 8) +
  labs(title = "KMeans Clustering (Mass vs Radius)")

### DBSCAN : 
# Load necessary libraries
library(dbscan)
ps_na = na.omit(selected_vars)


# 1. Scale the data (DBSCAN is distance-based, so scaling is essential)
scaled_data <- scale(ps_na)

# 2. Determine appropriate eps using kNN distance plot
kNNdistplot(scaled_data, k = 8)  # try k = minPts (usually 4–10)
abline(h = 1.5, col = "red", lty = 2)  # adjust threshold based on the elbow

# 3. Run DBSCAN
db <- dbscan(scaled_data, eps = 0.7, minPts = 8)  # adjust eps based on plot

# 4. Add cluster labels to original data
ps_na$cluster <- as.factor(db$cluster)  # 0 means noise

# 5. Visualize with PCA
fviz_cluster(list(data = scaled_data, cluster = db$cluster),
             geom = "point", ellipse = FALSE,
             main = "DBSCAN Clustering")
#### PCA : 





ps_pca_data = exo_koi_overlap %>% select(where(is.numeric))
dim(ps_pca_data)
extended_features <- c(
  "pl_rade", "pl_bmasse", "pl_orbper", "pl_orbeccen", "pl_insol", "pl_eqt",
  "st_teff", "st_rad", "st_mass", "st_logg", "st_met", "sy_dist",
  "pl_orbsmax", "st_metratio", "sy_vmag", "pl_eqterr1"  # example extras
)
selected_features <- c(
  "pl_rade",        # Planet radius (Earth)
  "pl_bmasse",      # Planet mass (Earth)
  "pl_orbper",      # Orbital period
  "pl_orbeccen",    # Eccentricity
  "pl_insol",       # Insolation flux
  "pl_eqt",         # Equilibrium temperature
  "st_teff",        # Stellar effective temperature
  "st_rad",         # Stellar radius
  "st_mass",        # Stellar mass
  "st_logg",        # Stellar surface gravity
  "st_met",         # Stellar metallicity
  "sy_dist"         # Distance
)
# Subset and remove rows with NA
ps_pca_data <- ps_pca_data[, selected_features]
ps_pca_data <- na.omit(ps_pca_data)
dim(ps_pca_data)


library(plotly)
library(factoextra)
head(ps_na)
# Assuming ps_na is already cleaned and numeric
# 1. Scale the data
scaled_vars <- scale(ps_pca_data)

# 2. Run PCA
pca_result <- prcomp(scaled_vars, center = T, scale. = T)

# 3. Get first 3 PCs
pc_data <- as.data.frame(pca_result$x[, 1:3])
colnames(pc_data) <- c("PC1", "PC2", "PC3")

# (Optional) Add labels or other grouping variables if you have them
# pc_data$cluster <- ps_na$cluster  # if you already have clusters

# 4. 3D plot using plotly
fig <- plot_ly(pc_data, x = ~PC1, y = ~PC2, z = ~PC3,
               type = 'scatter3d', mode = 'markers',
               marker = list(size = 3, color = ~PC1, colorscale = 'Viridis'))

fig <- fig %>% layout(title = "PCA: First 3 Principal Components",
                      scene = list(
                        xaxis = list(title = "PC1"),
                        yaxis = list(title = "PC2"),
                        zaxis = list(title = "PC3")
                      ))

fig


ps_pca_data$planet_type <- with(ps_pca_data, 
                                ifelse(pl_bmasse < 2 & pl_rade < 1.5, "Earth-like",
                                       ifelse(pl_bmasse >= 2 & pl_bmasse < 10 & pl_rade >= 1.5 & pl_rade < 2.5, "Super-Earth",
                                              ifelse(pl_bmasse >= 10 & pl_bmasse < 50 & pl_rade >= 2.5 & pl_rade < 4, "Neptune-like",
                                                     ifelse(pl_bmasse >= 50 & pl_rade >= 4, "Jupiter-like", 
                                                            "Other")))))
table(ps_pca_data$planet_type)


library(ggplot2)

ggplot(ps_pca_data, aes(x = pl_rade, y = pl_bmasse, color = planet_type)) +
  geom_point(alpha = 0.7) +
  scale_x_log10() + scale_y_log10() +
  labs(x = "Planet Radius (Earth Radius)", y = "Planet Mass (Earth Mass)", color = "Planet Type") +
  theme_minimal()

library(dbscan)

# Prepare the data matrix with only radius and mass
data_dbscan <- ps_pca_data[, c("pl_rade", "pl_bmasse")]

# Optional: scale features (recommended for DBSCAN)
data_scaled <- scale(data_dbscan)
# 2. Determine appropriate eps using kNN distance plot
kNNdistplot(data_scaled, k = 5)  # try k = minPts (usually 4–10)
abline(h = 0.4, col = "red", lty = 2)  # adjust threshold based on the elbow

# Run DBSCAN (you need to tune eps and minPts parameters)
set.seed(123)
dbscan_result <- dbscan(data_scaled, eps = 0.5, minPts = 7)

# Add cluster labels to your data frame
ps_pca_data$dbscan_cluster <- as.factor(dbscan_result$cluster)

# View cluster counts
table(ps_pca_data$dbscan_cluster)


### KMEANS : 
# Select the two features
data_kmeans <- ps_pca_data[, c("pl_rade", "pl_bmasse")]

# Remove rows with NA (if any)
data_kmeans <- na.omit(data_kmeans)

# Scale the data (recommended for KMeans)
data_scaled <- scale(data_kmeans)

# Set the number of clusters (e.g., 4)
set.seed(123)
kmeans_result <- kmeans(data_scaled, centers = 4, nstart = 25)

# Add cluster labels to the original dataset (only the rows used in clustering)
ps_pca_data$kmeans_cluster <- NA
ps_pca_data$kmeans_cluster[!is.na(ps_pca_data$pl_rade) & !is.na(ps_pca_data$pl_bmasse)] <- as.factor(kmeans_result$cluster)

# View cluster sizes
table(ps_pca_data$kmeans_cluster)


library(ggplot2)

ps_pca_data$kmeans_cluster <- as.factor(ps_pca_data$kmeans_cluster)

ggplot(ps_pca_data, aes(x = pl_rade, y = pl_bmasse, color = kmeans_cluster)) +
  geom_point(alpha = 0.7) +
  scale_x_log10() + scale_y_log10() +
  scale_color_manual(values = c("1" = "blue", "2" = "red", "3" = "darkgreen", "4" = "pink")) +
  labs(x = "Planet Radius (Earth Radius)", y = "Planet Mass (Earth Mass)", color = "Cluster") +
  theme_minimal()
#### Hierarchical clustering : 


# Assuming ps_pca_data is your data frame with numeric features pl_rade and pl_bmasse
data_for_clust <- ps_pca_data[, c("pl_rade", "pl_bmasse")]

# Scale the data (important for clustering)
scaled_data <- scale(data_for_clust)

# Compute distance matrix
dist_matrix <- dist(scaled_data)

# List of linkage methods to try
linkages <- c("complete", "average", "single", "ward.D2")

# Plot dendrograms for each linkage method
par(mfrow = c(2, 2))  # 2x2 plot layout

for(linkage in linkages) {
  hc <- hclust(dist_matrix, method = linkage)
  plot(hc, main = paste("Hierarchical Clustering -", linkage, "linkage"),
       xlab = "", sub = "", cex = 0.6)
}


library(ggplot2)

data_for_clust <- ps_pca_data[, c("pl_rade", "pl_bmasse")]
scaled_data <- scale(data_for_clust)
dist_matrix <- dist(scaled_data)

linkages <- c("complete", "average", "single", "ward.D2")

par(mfrow = c(2, 2))  # Optional: to visualize dendrograms if needed

for(linkage in linkages) {
  hc <- hclust(dist_matrix, method = linkage)
  clusters <- cutree(hc, k = 4)
  
  # Add cluster info to data frame
  plot_data <- data.frame(data_for_clust, cluster = as.factor(clusters))
  
  # Plot clusters
  p <- ggplot(plot_data, aes(x = pl_rade, y = pl_bmasse, color = cluster)) +
    geom_point(alpha = 0.7) +
    scale_x_log10() + scale_y_log10() +
    labs(title = paste("Hierarchical Clustering -", linkage, "linkage"),
         x = "Planet Radius (Earth Radius)", y = "Planet Mass (Earth Mass)",
         color = "Cluster") +
    theme_minimal() +
    scale_color_brewer(palette = "Set1")
  
  print(p)
}

#############
library(ggplot2)
library(RColorBrewer)

# Select the variables for clustering
selected_features <- c(
  "pl_rade", "pl_bmasse", "pl_orbper", "pl_orbeccen", "pl_insol", "pl_eqt",
  "st_teff", "st_rad", "st_mass", "st_logg", "st_met", "sy_dist"
)

# Extract and scale data
clust_data <- ps_pca_data[, selected_features]
scaled_data <- scale(clust_data)

# Compute distance matrix
dist_matrix <- dist(scaled_data)

# Linkage methods to trys
linkages <- c("complete", "average", "single", "ward.D2")

# Loop over linkage methods and plot clusters
for(linkage in linkages) {
  hc <- hclust(dist_matrix, method = linkage)
  clusters <- cutree(hc, k = 4)
  
  plot_data <- data.frame(
    pl_rade = ps_pca_data$pl_rade,
    pl_bmasse = ps_pca_data$pl_bmasse,
    cluster = as.factor(clusters)
  )
  
  p <- ggplot(plot_data, aes(x = pl_rade, y = pl_bmasse, color = cluster)) +
    geom_point(alpha = 0.7) +
    scale_x_log10() +
    scale_y_log10() +
    labs(
      title = paste("Hierarchical Clustering with", linkage, "linkage"),
      x = "Planet Radius (Earth Radius)",
      y = "Planet Mass (Earth Mass)",
      color = "Cluster"
    ) +
    theme_minimal() +
    scale_color_brewer(palette = "Set1")
  
  print(p)
}
#### more features : 
# Extended list of features (add as many numeric/continuous features as you want)
extended_features <- c(
  "pl_rade", "pl_bmasse", "pl_orbper", "pl_orbeccen", "pl_insol", "pl_eqt",
  "st_teff", "st_rad", "st_mass", "st_logg", "st_met", "sy_dist",
  "pl_orbsmax", "st_metratio", "sy_vmag", "pl_eqterr1"  # example extras
)

# Extract and scale
clust_data <- ps_pca_data[, extended_features]
clust_data <- na.omit(clust_data)  # remove rows with NA
scaled_data <- scale(clust_data)

# Distance matrix
dist_matrix <- dist(scaled_data)

# Linkage methods (you can pick one or try multiple as before)
linkage <- "ward.D2"
hc <- hclust(dist_matrix, method = linkage)
clusters <- cutree(hc, k = 4)

# Plotting clusters only on radius vs mass for visualization
plot_data <- data.frame(
  pl_rade = ps_pca_data$pl_rade[rownames(clust_data)],
  pl_bmasse = ps_pca_data$pl_bmasse[rownames(clust_data)],
  cluster = as.factor(clusters)
)

library(ggplot2)
ggplot(plot_data, aes(x = pl_rade, y = pl_bmasse, color = cluster)) +
  geom_point(alpha = 0.7) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    title = paste("Hierarchical Clustering with", linkage, "linkage (More features)"),
    x = "Planet Radius (Earth Radius)",
    y = "Planet Mass (Earth Mass)",
    color = "Cluster"
  ) +
  theme_minimal() +
  scale_color_brewer(palette = "Set1")

