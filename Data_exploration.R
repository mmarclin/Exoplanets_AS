### KOI_table

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
summary(KOI_table)
pairs(KOI_table)
attach(KOI_table)






detach(KOI_table)


