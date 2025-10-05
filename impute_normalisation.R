library(missForest)
library(openxlsx)

# Example dataset with missing values
data(iris)
set.seed(111)
iris.na <- iris
for (i in 1:4) iris.na[sample(150, sample(20, 1)), i] <- NA

# Impute missing values using missForest
iris.imputed <- missForest(iris.na)$ximp

# View the imputed dataset
head(iris.imputed)

dtIN <- read.xlsx("proteinGroups_intensity_p-Value-posthoc-annotation.xlsx")

head(dtIN)
dt_intensity <- dtIN[, grep("Intensity", colnames(dtIN))]
dt_intensity[dt_intensity == 0] <- NA

dt_intensity.imputed <- missForest(dt_intensity)$ximp

saveRDS(dt_intensity.imputed, "/Users/163726/Library/CloudStorage/OneDrive-Personal/SideHustle/Mai/Aurantiamide/randomForest_imputed_intensity.rds")

head(dt_intensity.imputed)
head(dt_intensity[1:100,])


