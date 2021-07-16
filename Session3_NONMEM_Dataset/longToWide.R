# Convert long data simulated from Session 2 to wide data
# Generate PK only and PKPD dataset for demonstration

library(tidyverse)
library(data.table)

# read in data from Session2_Simulations
PKPD <- read_csv("Session2_Simulations/r4pmx_dataset_simulation.csv", na=".")

# convert long to wide data
PKPD_wide <- PKPD %>% 
  mutate(type = if_else(type==1, "PK", "PD")) %>%   
  select(-tad, -cmt, -evid, -mdv, -dose, -bmi) %>%
  pivot_wider(names_from = c(type, nom_tad), values_from = dv) # convert long to wide

PKPD_wide <- PKPD_wide %>% 
  group_by(ID) %>% 
  fill(names(PKPD_wide), .direction = "downup") %>% 
  filter(row_number() == 1) %>% 
  select(-PK_0)

# rename covariates
PKPD_wide <- PKPD_wide %>% 
  mutate(sex = if_else(sex == 1, "Female", "Male"),
         ethnicity = if_else(ethnicity == 1, "Hispanic", "Non-Hispanic"),
         geno = case_when(geno == 1 ~ "wild", geno == 2 ~ "heterozygous", geno == 3 ~ "homozygous"))

PK_wide <- PKPD_wide[,c(1:11,grep("PK",names(PKPD_wide)))]

write.csv(PK_wide, "data/PK_Data_Session3.csv", row.names = F)
write.csv(PKPD_wide, "data/PKPD_Data_Session3.csv", row.names = F)