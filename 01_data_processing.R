library(missForest)
library(openxlsx)

# Example dataset with missing values
data(iris)
set.seed(111)
iris.na <- iris
for (i in 1:4) iris.na[sample(150, sample(20, 1)), i] <- NA


dtIN <- read.xlsx("proteinGroups_intensity_p-Value-posthoc-annotation.xlsx")

head(dtIN)
dt_intensity <- dtIN[, grep("Intensity", colnames(dtIN))]
dt_intensity[dt_intensity == 0] <- NA



