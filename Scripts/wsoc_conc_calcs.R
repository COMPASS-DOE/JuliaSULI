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


## 4. Calculate Kit Replicate % Error and Average Concentrations ---------------

# make list of replicate kits
rep_list <- wsom_extracts %>%
  filter(grepl("_rep", sample_name)) %>%
  select(sample_name)

# find the other data points for the replicate sites and bind together
kit_replicates <- wsom_extracts %>%
  filter(substr(wsom_extracts$sample_name,1,10) %in% 
                substr(rep_list$sample_name,1,10)) %>%
  select(kit_id, wsoc_conc, wson_conc, transect_location)

# average calcs
kit_replicates_averages <- kit_replicates %>%
  mutate(kit_id = substr(kit_id, 1,4)) %>%# replicate kits have same name
  group_by(kit_id, transect_location) %>%
  summarize(avg_wsoc_conc = mean(wsoc_conc),
            avg_wson_conc = mean(wson_conc))

# calculate percent difference
npoc_tdn_perc_diff <- kit_replicates_averages %>%
  select(kit_id, transect_location)

npoc_tdn_perc_diff$wsoc_perc_error <- NA # create columns for calcs
npoc_tdn_perc_diff$wson_perc_error <- NA

for (i in 1:nrow(rep_list)) {
  
  # filter the replicate dataset to just the two observations that match
  kit_replicates_onekit <- kit_replicates %>%
    filter(substr(kit_id,1,4) == npoc_tdn_perc_diff$kit_id[i])
  
  npoc_tdn_perc_diff$wsoc_perc_error[i] <-     # (second - initial)/initial
    abs(((kit_replicates_onekit[grepl("_rep", kit_replicates_onekit$kit_id),
                          "wsoc_conc"] - 
      kit_replicates_onekit[!grepl("_rep", kit_replicates_onekit$kit_id), 
                             "wsoc_conc"]) /
      kit_replicates_onekit[!grepl("_rep", kit_replicates_onekit$kit_id), 
                             "wsoc_conc"]) * 100)
  
  npoc_tdn_perc_diff$wson_perc_error[i] <-     
    abs(((kit_replicates_onekit[grepl("_rep", kit_replicates_onekit$kit_id),
                          "wson_conc"] - 
      kit_replicates_onekit[!grepl("_rep", kit_replicates_onekit$kit_id), 
                            "wson_conc"]) /
      kit_replicates_onekit[!grepl("_rep", kit_replicates_onekit$kit_id), 
                            "wson_conc"]) * 100)
}

# write out the perc error calcs
write.csv(npoc_tdn_perc_diff, "./wsom_extracts_npoc_tdn_perc_error_calcs.csv")
