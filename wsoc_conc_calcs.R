## This script reads in fully processed NPOC/TDN data and a csv file with
## WSOM soil extraction procedure notes to calculate WSOC and WSON 
## concentrations normalized to the amount of water and soil used in the
## extraction process.
##
## Data are saved in the JuliaSULI repository.

## Created: 2022/11/22
## Julia McElhinny

## 1. Load Packages and Data ---------------------------------------------------

require(tidyverse)

# check to make sure working directory is main repository folder
getwd()

# load in csv files
wsoc_masses <- read.csv("./wsoc_extracts_masses.csv", na.strings = "N/A")
wsoc_npoc_tdn <- read.csv("./EC1_WSOC_Extracts_NPOC_TDN_L0B_20221118.csv")


## 2. Merge Data Files ---------------------------------------------------------

# make sure replicate kits are named that way
for (n in 1:nrow(wsoc_npoc_tdn)) {
  old_name <- wsoc_npoc_tdn[n,"kit_id"]
  
  # find the rows with "_rep" in the sample name
  if (grepl("_rep", wsoc_npoc_tdn[n,1])) {
    
    wsoc_npoc_tdn[n,"kit_id"] <- paste0(old_name, "_rep")
  }
}

# merge by kit_id
wsom_extracts_withblanks <- left_join(wsoc_npoc_tdn, wsoc_masses, by = 
                                        c("kit_id", "transect_location"))

# remove blanks (npoc/tdn data has already been blank corrected)
wsom_extracts <- wsom_extracts_withblanks %>%
  filter(!grepl("blank", sample_name))


## 3. Calculate Conc/g Dry Weight ----------------------------------------------
## formula follows the normalization methods of Wardinski et al. 2022
## (mg/L NPOC * L water used) / g soil used

# set value of water density
water_density <- 0.99823 # density of water @ room temp (g/mL)
water_density_gL <- water_density*1000

# calculate L of water used for extraction
wsom_extracts$water_L <- wsom_extracts$water_mass / water_density_gL

# apply formula to get mg WSOM/g soil of carbon and nitrogen
wsom_extracts <- wsom_extracts %>%
  mutate(wsoc_conc = (npoc_mgl*water_L)/soil_mass,
         wson_conc = (tdn_mgl*water_L)/soil_mass)


## 4. Calculate Percent Error with Replicates ----------------------------------

# make list of replicate kits
rep_list <- wsom_extracts %>%
  filter(grepl("_rep", kit_id)) %>%
  select(kit_id, wsoc_conc, wson_conc)
