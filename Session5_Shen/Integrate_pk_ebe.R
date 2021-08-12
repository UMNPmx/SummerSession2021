#library
library(tidyverse)
library(patchwork)

#setwd
setwd(".")
getwd()

#read dataset
df <- read_csv("r4pmx_dataset_simulation1.csv")
dat <- read.table (file='sdtab003', skip=1, header=T)

#param
param <- dat %>% select(ID, KTR, CL, VC) %>% distinct(ID, .keep_all = TRUE)

#merge param
df_new <- merge(df, param, by="ID", all=TRUE)

#write.csv
df_new <- df_new %>% select(C, everything()) %>% mutate(C = ".")

write.csv(df_new, "r4pmx_dataset_simulation1_pkparam.csv", row.names = FALSE)
