setwd("/home/midus_1114420/edwa0506")

library(tidyverse)

#######################
# Explore MIDJA files
########################

MIDJA_meta <- read_csv("data/Imputed_SNPs/MIDJA Phenotypic Data/MIDJAGeneticsMetadata_n=328_07-12-22.csv")

head(MIDJA_meta)

MIDJA_meta  |>
  select(-M2MRMJ_ID_S) |>      
  map(unique)
# MIDJA is all Japanes from MIDJA2. All use Omni10 array
# Other MIDJA files just include phenotyes

#######################
# Explore M2 files
########################

M2_meta <-read_csv("data/Imputed_SNPs/MIDUS 2 Metadata/M2GeneticsMetadata_n=980_07-12-22.csv")

head(M2_meta)

M2_meta  |>
  select(-M2MRMJ_ID_S) |>      
  map(unique)
# M2 includes a range of self-reported races and omni10 and omni11

nrow(M2_meta)
# 980 people

M2_meta  |>
  select(RACE, OMNICHIP) |>
  table() 
# 1 refused. 1 Don't know. 22 Other.
# 772 Whites
# of Whites 137 Omni11 and 635 Omni10 

###########################
# Explore M Refresher files
###########################

MR_meta <-read_csv("data/Imputed_SNPs/MIDUS Refresher Metadata/MRGeneticsMetadata_n=809_07-12-22.csv")

head(MR_meta)
nrow(MR_meta)
#809

MR_meta  |>
  select(-M2MRMJ_ID_S) |>      
  map(unique)
# Only omni11

MR_meta  |>
  select(RACE) |>
  table() 
#575 White
# 3 Don't know
# 50 other
# 2 refused 


