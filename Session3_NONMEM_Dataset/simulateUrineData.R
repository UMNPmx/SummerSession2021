library(mrgsolve)
library(truncnorm) #functionality related to truncated nomral distribution
library(tidyverse)
theme_set(theme_bw())

# Mod File
pkpdurine_model <- '

[set] delta = 0.01, end = 50 

[PROB]
 - author: Shen Cheng, Jay Wen
 - date: `r Sys.Date()`
 - comments: this is a demonstration
 - PK: one compartment model + transit compartment absportion + first order elimination
 - PD: indirect response model with drug inhibiting `kin`
 - Urine: collect urine output
 - allometric scaling on flow and volume based parameters
 - random effects added on `ktr`, `clr`, `clnr`, `vc`, `ic50` and `kin` as exponential distributions
 - proportional error model
 - model covariates: geno ~ clnr, wt and age ~ clr and sex ~ kin

[PARAM] @annotated
tvclr     :  0.7 : typical value of renal clearance (L/hr)
tvclnr    :  0.3 : typical value of non-renal clearance (L/hr)
tvvc      :  20  : typical value of central volume of distribution (L)
tvktr     :  1   : typical value of transit rate (/hr) 
tvic50    :  1   : typical value of concentration which introduce 50% of inhibitory effects (mg/L)
tvkin     :  50  : synthetic rate (/hr)
tvkout    :  5   : degradation rate (/hr)
geno2clnr :  0.5 : genotype 2 effect on clnr (IM)
geno3clnr :  0.2 : genotype 3 effect on clnr (PM)
ageclr    :  -0.3: age effect on clr
sexkin    :  0.8 : sex effect on kin
geno      :  1   : reference genotype (extensive metabolizor: 1/intermediate metabolizor: 2/poor metabolizor: 3)
age       :  50  : reference age (years)
sex       :  0   : reference gender (female: 0/male: 1)
wt        :  70  : refernce wt (kg)

[CMT]  @annotated 
depot  : depot compartment (mass)
trans1 : transit compartment 1 (mass)
trans2 : transit compartment 2 (mass)
trans3 : transit compartment 3 (mass)
cent   : central compartment (mass)
resp   : response compartment (conc)
urine  : urine compartment (mass)

[MAIN]
double clr = tvclr*pow(wt/70, 0.75)*pow(age/50, ageclr)*exp(ETA(1));
double clnr = tvclnr*exp(ETA(2));
  if (geno == 2) clnr = tvclnr*exp(ETA(2))*geno2clnr;
  if (geno == 3) clnr = tvclnr*exp(ETA(2))*geno3clnr;
double vc = tvvc*pow(wt/70,1)*exp(ETA(3)); 
double ktr = tvktr*exp(ETA(4));
double ic50 = tvic50*exp(ETA(5));
double kin = tvkin*pow(sexkin, sex)*exp(ETA(6));
double kout = tvkout; 

resp_0 = kin/kout; //initialize the response compartment at steady state

[OMEGA] @annotated @correlation
iivclr  : 0.2 : interindividual variability (iiv) on clr
iivclnr : 0.1 0.2 : interindividual variability (iiv) on clnr (10% correlation between clr and clnr)
iivvc   : 0.25 0.25  0.2 : iiv on vc (25% correlation between clr/clnr and vc)

[OMEGA] @annotated
iivktr  : 0.1 : iiv on ktr
iivic50 : 0.3 : iiv on ic50
iivkin  : 0.2 : iiv on kin

[SIGMA] @annotated
ruvpk   : 0.01 : ruv on pk obs
ruvpd   : 0.04 : ruv on pd obs

$ODE
//PK 
double cp   = cent/vc; 

dxdt_depot  = -ktr*depot;
dxdt_trans1 = ktr*(depot - trans1);
dxdt_trans2 = ktr*(trans1 - trans2);
dxdt_trans3 = ktr*(trans2 - trans3);
dxdt_cent   = ktr*trans3 - (clr+clnr)*cp;
dxdt_urine  = clr*cp;

//PD
double inh = cp /(cp + ic50);
double res = resp;

dxdt_resp = kin*(1-inh) - kout*resp; 

[TABLE]
capture pkobs = cp*(1 + EPS(1));
capture urineobs = urine*(1 + EPS(1)); 
capture pdobs = resp*(1 + EPS(2)); 

'

# compile code
pkpdurine <- mcode("pkpdurine_model", pkpdurine_model)

# Simulate dataset
set.seed(0712) #set a seed for reproducibility 
data <- 
  data.frame(ID = 1:100,                                              #100 subjects
             amt = sample(c(50,100,200),replace = TRUE, size=100),    #sample 100 doses from c(50, 100, 200) with replacement
             cmt = 1,                                                 #dose in compartment 1
             evid = 1,                                                #event id 1: dose event
             time = 5,                                                #non-zero dose time to accommodate predose sample
             geno = sample(c(1,2,3), replace = TRUE, size=100),       #genotype (1: EM/2: IM/3: PM)
             age = rtruncnorm(n=100, a=0, b=Inf, mean=50, sd=5),      #age (trunc N(50,5) (lower bound 0)
             sex = sample(c(0,1), replace = TRUE, size=100),          #sex (0: female/1: male)
             ethnicity = sample(c(0,1), replace = TRUE, size=100)     #ethnicity (0: hispanic/1: non-hispanic)
  )

#generate weight (wt), height (ht) and creatine clearance (crcl) 
#assuming these parameters are correlated with sex

output <- data.frame(NULL)
for (i in 1:100){
  temp <- subset(data, ID == i)  #subset one individual for one time
  temp <- 
    temp %>% 
    mutate(wt = if_else(sex == 0,
                        rtruncnorm(n=1, a=0, b=Inf, mean=60, sd=10), #wt female (trunc N(60,10) (lower bound 0) (kg)
                        rtruncnorm(n=1, a=0, b=Inf, mean=70, sd=10)) #wt male (trunc N(70,10) (lower bound 0) (kg)
    ) %>%
    mutate(ht = if_else(sex == 0,
                        rtruncnorm(n=1, a=0, b=Inf, mean=160, sd=10), #ht female (trunc N(160,10) (lower bound 0) (cm)
                        rtruncnorm(n=1, a=0, b=Inf, mean=175, sd=10)) #ht male (trunc N(175,10) (lower bound 0) (cm)
    ) %>% 
    
    mutate(crcl = if_else(sex == 0,
                          rtruncnorm(n=1, a=0, b=Inf, mean=115, sd=10), #crcl female (trunc N(115,10) (lower bound 0) (ml/min)
                          rtruncnorm(n=1, a=0, b=Inf, mean=130, sd=10)) #crcl male (trunc N(130,10) (lower bound 0) (ml/min)
    )
  output <- rbind(output, temp) #row bind temp dataset after processing
}

# Simulation
# concentration
set.seed(0712)
sim <- 
  pkpdurine %>% 
  data_set(data) %>% 
  mrgsim(end=100, delta=0.5, tad=TRUE, obsonly=TRUE) %>% #`tad=TRUE` allows the calculation of the "time after dose"
  as.data.frame()

# Capture observation data
obs <- 
  sim %>% 
  filter(tad %in% c(-2, 0.5, 1, 2, 4, 8, 12, 24, 48)) %>% #select datapoints based on study design from simulated dataset
  select(ID, time, tad, pkobs, urineobs, pdobs) #only keep these five columns

## visual check
ggplot(obs, aes(x=tad, y=pkobs, group=ID))+
  geom_point() + geom_line() +
  ylab("Plasma concentration")

ggplot(obs, aes(x=tad, y=urineobs, group=ID))+
  geom_point() + geom_line() +
  ylab("Cumulative amount in urine")

# urine is cumulative, make it amount in each interval: 8, 24, 48
urineObs <- obs %>% 
  filter(tad %in% c(8, 24, 48)) %>%
  group_by(ID) %>% 
  mutate(urinediff = if_else(is.na(urineobs  - lag(urineobs)), urineobs, urineobs  - lag(urineobs))) %>% # calculate the difference
  select(ID, tad, urineobs, urinediff)

# combine pkobs, urineobs, pdobs
urineObsClean <- urineObs %>% 
  mutate(dv = urinediff, type = "urine") %>% 
  select(-urineobs, -urinediff)

pkObs <- obs %>% 
  mutate(dv = pkobs, type = "pk") %>% 
  select(ID, tad, dv, type)

pdObs <- obs %>% 
  mutate(dv = pdobs, type = "pd") %>% 
  select(ID, tad, dv, type)

obsAll <- bind_rows(pkObs, urineObsClean) %>% bind_rows(., pdObs) %>% arrange(ID)

# convert long to wide
obsAllWide <- obsAll %>% 
  pivot_wider(names_from = c(type, tad), 
              values_from = dv)

obsAllWide <- obsAllWide %>% 
  group_by(ID) %>% 
  fill(names(obsAllWide), .direction = "downup") %>% 
  slice(1)

# combine covariates with observations
PKPDUrineData <- full_join(obsAllWide, data[c("ID", "amt","geno", "age", "sex", "ethnicity")], by = "ID") %>% 
  select(-"pk_-2")

write_csv(PKPDUrineData, "data/PKPDUrine_Data_Session3.csv")