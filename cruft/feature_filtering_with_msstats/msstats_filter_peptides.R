if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("MSstats")

library('MSstats')
library(MSstatsConvert)
library(data.table)

rawinput <- read.csv(file="peptided_msstats.csv")
setDT(input)  
input = MSstatsPrepareForDataProcess(input, 2, NULL)
input_informative = MSstatsSelectFeatures(input, "topN") # feature selection

filtereddf = input_informative[input_informative$remove == 0, ] #look at peptides that weren't filtered out
View(filtereddf)

#protein ID from the paper:
# PGM1: P36871
#  RRS1: Q15050
# ACO2: Q99798
