### Install and load necessary libraries

# Load libraries ----
library(dplyr)       # Data manipulation
library(GGally)      # Pair plots
library(corrplot)    # Correlation heatmaps
library(ggplot2)     # Visualization
library(rgl)         # 3D visualization
library(mice)        # Handling missing data
library(VIM)         # Visualizing missing data
#####
### Import the data 
KOI_table_full = read.csv("C:/Users/antoi/OneDrive/Bureau/AppliedStats/PROJECT EXOPLANETS/KOI_table_2025.03.23_01.50.59.csv",comment.char = "#",header = TRUE, stringsAsFactors = T)

### Data Exploration
head(KOI_table_full)       # First few rows
names(KOI_table_full)      # Column names
str(KOI_table_full)        # Structure of the dataset
dim(KOI_table_full)        # Dataset dimensions
attach(KOI_table_full)



data = KOI_table_full
### Data Cleaning and Preprocessing

# First we need to create the dataset with the target variable and the features 
# the features are all the variables - those who are not useful for the model

columns_to_remove <- c("kepid", "kepoi_name", "kepler_name",
                       "koi_pdisposition","koi_score" ,
                       "koi_teq_err1", "koi_teq_err2", "koi_tce_delivname")

KOI_table = KOI_table_full %>% select(-all_of(columns_to_remove))

# Remove unconfirmed candidate exoplanets
KOI_table <- KOI_table %>% filter(koi_disposition != "CANDIDATE")

# Encode 'koi_disposition' as binary: 1 (Confirmed), 0 (False Positive)
KOI_table$koi_disposition <- ifelse(KOI_table$koi_disposition == "CONFIRMED", 1, 0)
# encode them as factors ?

# Handle missing data : 
# replace by mean / closest neighbors value ... 
na_par_var = numeric(ncol(KOI_table))
for(i in 1:ncol(KOI_table)){
  na_par_var[i] = sum(is.na(KOI_table[,i]))
}

na_par_var
Var_NA = colnames(KOI_table)[na_par_var!=0]
Var_NA

# 1st approach : remove rows with NA

KOI_table_clean <- na.omit(KOI_table)  # Remove rows with NA values

### Data Analysis : 

# Correlation matrix
cor_matrix <- cor(KOI_table_clean[, -which(names(KOI_table_clean) == "koi_disposition")], use = "complete.obs")
corrplot(cor_matrix, method = "color", tl.cex = 0.6, tl.srt = 45)
# some features are correlated maybe summarize some of them, PCA... 
# not too much correlation 

# Pairplot for selected features
selected_features <- KOI_table_clean[, c("koi_fpflag_nt", "koi_period", "koi_duration", "koi_teq", "koi_steff", "koi_disposition")]
ggpairs(selected_features, aes(color = as.factor(koi_disposition), alpha = 0.5))


# Boxplots (Standardized Data)
KOI_table_scaled <- KOI_table_clean
KOI_table_scaled[, -which(names(KOI_table_scaled) == "koi_disposition")] <- scale(KOI_table_scaled[, -which(names(KOI_table_scaled) == "koi_disposition")])
boxplot(KOI_table_scaled, las = 2, main = "Standardized Features Boxplot")

KOI_table_scaled

# Class distribution ----
barplot(prop.table(table(KOI_table_clean$koi_disposition)), col = c("red", "blue"), ylim = c(0,1), main = "Class Distribution", ylab = "Proportion", xlab = "KOI Disposition")
# Plot density of each features for each group : 
list_features = setdiff(colnames(KOI_table_clean),c("koi_disposition"))

for (var in list_features) {
  print(
    ggplot(KOI_table_clean, aes(x = !!sym(var), fill = factor(koi_disposition))) +
      geom_density(alpha = 0.5) +
      scale_fill_manual(values = c("steelblue4", "indianred", "darkgreen")) +
      theme_minimal() +
      ggtitle(paste("Density Plot of", var))
  )
}

# Our quantitative variables don't seem to be discriminating individually the classes too much 
var = koi_steff
ggplot(KOI_table_clean) +
  aes(x = koi_steff, fill = factor(koi_disposition))  +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("steelblue4", "indianred")) +
  theme_minimal()
# Density for 2 features at the same time: 
# Choose two numeric features
feature_x <- "koi_steff"
feature_y <- "koi_period"

ggplot(KOI_table_clean, aes_string(x = feature_x, y = feature_y, color = "koi_disposition")) +
  geom_density_2d() +
  theme_minimal() +
  labs(title = paste("2D Density Plot of", feature_x, "vs", feature_y),
       x = feature_x, y = feature_y)


# Loop over each class and plot 3D density
unique_classes <- unique(KOI_table_clean$koi_disposition)

for (cls in unique_classes) {
  # Subset and clean data
  data_subset <- KOI_table_clean %>%
    filter(koi_disposition == cls) %>%
    filter(!is.na(.data[[feature_x]]), !is.na(.data[[feature_y]]))
  
  # Perform 2D kernel density estimate
  kde <- kde2d(data_subset[[feature_x]], data_subset[[feature_y]], n = 100)
  
  # Create 3D surface plot
  plot <- plot_ly(x = kde$x, y = kde$y, z = kde$z, type = "surface") %>%
    layout(
      title = paste("3D Density for Class:", cls),
      scene = list(
        xaxis = list(title = feature_x),
        yaxis = list(title = feature_y),
        zaxis = list(title = "Density")
      )
    )
  
  print(plot)
}

# Now we can represent 2 on the same graphs: 

# Choisir les deux features
feature_x <- "koi_steff"
feature_y <- "koi_duration"

# Choisir deux classes à comparer (tu peux les modifier)
classes_to_plot <- unique(KOI_table_clean$koi_disposition)[1:2]

# Palette et opacités
colors <- c("red", "blue")
opacities <- c(1.0, 0.5)  # deuxième classe sera plus transparente

# Initialiser le plot
final_plot <- plot_ly()

# Boucle pour ajouter les surfaces
for (i in seq_along(classes_to_plot)) {
  cls <- classes_to_plot[i]
  
  # Sous-échantillon pour la classe
  data_subset <- KOI_table_clean %>%
    filter(koi_disposition == cls) %>%
    filter(!is.na(.data[[feature_x]]), !is.na(.data[[feature_y]]))
  
  # Calcul densité
  kde <- kde2d(data_subset[[feature_x]], data_subset[[feature_y]], n = 100)
  
  # Ajouter la trace
  final_plot <- final_plot %>% add_surface(
    x = kde$x,
    y = kde$y,
    z = kde$z,
    name = paste("Classe:", cls),
    showscale = FALSE,
    opacity = opacities[i],
    surfacecolor = matrix(1, nrow = length(kde$x), ncol = length(kde$y)),
    colorscale = list(c(0, 1), c(colors[i], colors[i]))
  )
}

# Ajouter le layout
final_plot <- final_plot %>%
  layout(
    title = paste("3D Densité conjointe de", feature_x, "et", feature_y),
    scene = list(
      xaxis = list(title = feature_x),
      yaxis = list(title = feature_y),
      zaxis = list(title = "Densité")
    ),
    legend = list(title = list(text = "Classes"))
  )

final_plot

#####


# We can also do some statistical tests on the mean to detect if some features 
# have different means depending on the classes : 

# Are our features normally distributed : 
#  perform Shapiro-Wilk test on  each numeric column 
normality_tests <- sapply(KOI_table_clean, function(x) {
  if (is.numeric(x)) {
    shapiro_result <- shapiro.test(x)
    return(shapiro_result$p.value)  # Extract p-value from the test result
  } else {
    return(NA)  # If the column is not numeric, return NA
  }
})
# We have too many data to perform this test :  we can try anderson's darling test  
install.packages("nortest")
library(nortest)
anderson_tests <- sapply(KOI_table_clean, function(x) {
  if (is.numeric(x)) {
    ad_result <- ad.test(x)  # Anderson-Darling test
    return(ad_result$p.value)
  } else {
    return(NA)
  }
})

# Print the p-values from Anderson-Darling test
anderson_tests

# Histogram : 
for (col in colnames(KOI_table_clean)) {
  if (is.numeric(KOI_table_clean[[col]])) {
    hist(KOI_table_clean[[col]], prob = TRUE, main = paste("Histogram of", col), xlab = col)
    curve(dnorm(x, mean = mean(KOI_table_clean[[col]]), sd = sd(KOI_table_clean[[col]])),
          add = TRUE, col = "blue", lwd = 2)
  }
}

# Q-Q plots for each numeric column
for (col in colnames(KOI_table_clean)) {
  if (is.numeric(KOI_table_clean[[col]])) {
    qqnorm(KOI_table_clean[[col]], main = paste("Q-Q plot of", col))
    qqline(KOI_table_clean[[col]], col = "red")
  }
}
# => some features don't follow a normal
# distribution : student t test for these variables are not appropriate 

## Test d'égalité des moyennes de chaque classe sur features : 
#Condition : Ce test suppose que les données suivent une distribution normale. 
#Si cette hypothèse est violée, vous devrez utiliser un test non paramétrique comme le test de Wilcoxon.
t.test(KOI_table_clean$koi_period~KOI_table_clean$koi_disposition)
t.test(KOI_table_clean$koi_steff~KOI_table_clean$koi_disposition)
t.test(KOI_table_clean$koi_impact~KOI_table_clean$koi_disposition)
t.test(KOI_table_clean$ra~KOI_table_clean$koi_disposition)

## Non parametric Wilcoxon Test : 
# for the equality of the median 
# H0 : median are equal in both groups 
wilcox.test(KOI_table_clean$koi_period~KOI_table_clean$koi_disposition)
# median are not equal 

### Analysis of the link between 2 qualitative variables : 
# Create a contingency table
t <- table(KOI_table_clean$koi_fpflag_ss, KOI_table_clean$koi_disposition)
t = t[,2:3]
prop_t <- prop.table(t, margin = 1) 
prop_t


# Calculate column proportions (proportions for each class of koi_disposition)
cprop_t <- prop.table(t, margin = 2)  # Proportion by column (koi_disposition)
cprop_t


# TEst chi-squared of independance between the 2 categorical variables
chi_test <- chisq.test(t)
# h0 : they are independant 
chi_test
# we can reject H0 => the 2 variables are not independant : 
# we can do this test for every qualitative variable and 
# remove the ones that are independant because they don't 
# bring informations about the target variable 



### PCA Analysis

# since we have a bit of correlation in our features we can do a PCA to reduce it  : 
pca <- prcomp(KOI_table_clean[, -which(names(KOI_table_clean) == "koi_disposition")], scale. = TRUE)
a = summary(pca)
eigenvalues <- pca$sdev^2
barplot(eigenvalues)
plot(cumsum(eigenvalues))

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



KOI_table_PCA = pca$x
# On ne garde que les p premières composantes : 
p = 20 
# Pour ..
cumsum(eigenvalues)[20]/tail(cumsum(eigenvalues),1)
# % de la variance totale
KOI_table_PCA = pca$x[,1:20]
dim(KOI_table_PCA)

### Models prediction




# we can do : logistic regression / SVMs / RandomForest 
### train models and do PCA with scaled data

