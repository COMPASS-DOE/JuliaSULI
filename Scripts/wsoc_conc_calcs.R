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
wsom_masses <- read_csv("./WSOM Extracts Data/wsoc_extracts_masses.csv")
wsom_npoc_tdn <- read_csv("./WSOM Extracts Data/EC1_WSOM_Extracts_NPOC_TDN_L0B.csv")

# read in soils RDS file for GWC
soils_all <- readRDS("./Downloaded Data/soils_data_merged_withmeta.rds")


## 2. Merge Data Files ---------------------------------------------------------

# pull out necessary GWC values for back calculation to field moist soils
extracts_gwc <- soils_all %>%
  filter(kit_id %in% substr(wsom_npoc_tdn$kit_id,1,4) & 
           transect_location %in% wsom_npoc_tdn$transect_location) %>%
  select(kit_id, transect_location, gwc_perc)

# add GWC to the npoc/tdn data frame
wsom_npoc_tdn <- left_join(wsom_npoc_tdn, extracts_gwc, by =
                             c("kit_id", "transect_location"))

# make sure replicate kits are named that way
for (n in 1:nrow(wsom_npoc_tdn)) {
  
  old_name <- wsom_npoc_tdn[n,"kit_id"]
  
  # find the rows with "_rep" in the sample name
  if (grepl("_rep", wsom_npoc_tdn[n,4])) {
    
    wsom_npoc_tdn[n,"kit_id"] <- paste0(old_name, "_rep")
  }
}

# merge by kit_id
wsom_extracts_withblanks <- left_join(wsom_npoc_tdn, wsom_masses, by = 
                                        c("kit_id", "transect_location"))

# remove blanks (npoc/tdn data has already been blank corrected)
wsom_extracts <- wsom_extracts_withblanks %>%
  filter(!grepl("blank", sample_name))


## 3. Calculate Conc/g Dry Weight ----------------------------------------------
## formula follows the normalization methods of Wardinski et al. 2022
## (mg/L NPOC * L water used) / g soil used

# set value of water density
water_density <- 0.99777 # density of water @ room temp 22C (g/mL)
water_density_gL <- water_density*1000

# calculate L of water used for extraction
wsom_extracts <- wsom_extracts %>%
  mutate(water_mass = as.numeric(water_mass),
         soil_mass = as.numeric(soil_mass)) %>%
  mutate(water_L = water_mass / water_density_gL)

# apply formula to get mg WSOM/g dry soil of carbon and nitrogen
wsom_extracts <- wsom_extracts %>%
  mutate(wsoc_conc = (npoc_mgl*water_L)/soil_mass,
         wson_conc = (tdn_mgl*water_L)/soil_mass)

## also calculate normalized to field moist soil (back calculate from GWC)
# calculate the field moist equivalent of the freeze dried soil used
wsom_extracts <- wsom_extracts %>%
  mutate(field_soil = ((gwc_perc/100)*soil_mass)+soil_mass)

# apply formula for mg WSOM/g field moist soil
wsom_extracts <- wsom_extracts %>%
  mutate(wsoc_conc_field = (npoc_mgl*water_L)/field_soil,
         wson_conc_field = (tdn_mgl*water_L)/field_soil)

# write out full concentrations data frame as csv
#write_csv(wsom_extracts, "./WSOM Extracts Data/wsom_extracts_wsoc_wson_conc.csv")


## 4. Calculate Kit Replicate % Error and Average Concentrations ---------------

# make list of replicate kits
rep_list <- wsom_extracts %>%
  filter(grepl("_rep", sample_name)) %>%
  select(sample_name)

# find the other data points for the replicate sites and bind together
kit_replicates <- wsom_extracts %>%
  filter(substr(wsom_extracts$sample_name,1,10) %in% 
                substr(rep_list$sample_name,1,10)) %>%
  select(kit_id, wsoc_conc, wson_conc, transect_location, wsoc_conc_field, 
         wson_conc_field)

# average calcs
kit_replicates_averages <- kit_replicates %>%
  mutate(kit_id = substr(kit_id, 1,4)) %>%# replicate kits have same name
  group_by(kit_id, transect_location) %>%
  summarize(avg_wsoc_conc = mean(wsoc_conc),
            avg_wson_conc = mean(wson_conc),
            avg_wsoc_conc_field = mean(wsoc_conc_field),
            avg_wson_conc_field = mean(wson_conc_field))

# calculate percent difference
npoc_tdn_perc_diff <- kit_replicates_averages %>%
  select(kit_id, transect_location)

npoc_tdn_perc_diff$wsoc_perc_error <- NA # create columns for calcs
npoc_tdn_perc_diff$wson_perc_error <- NA
npoc_tdn_perc_diff$wsoc_field_perc_error <- NA
npoc_tdn_perc_diff$wson_field_perc_error <- NA

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
  
  npoc_tdn_perc_diff$wsoc_field_perc_error[i] <-     
    abs(((kit_replicates_onekit[grepl("_rep", kit_replicates_onekit$kit_id),
                                "wsoc_conc_field"] - 
            kit_replicates_onekit[!grepl("_rep", kit_replicates_onekit$kit_id), 
                                  "wsoc_conc_field"]) /
           kit_replicates_onekit[!grepl("_rep", kit_replicates_onekit$kit_id), 
                                 "wsoc_conc_field"]) * 100)
  
  npoc_tdn_perc_diff$wson_field_perc_error[i] <-     
    abs(((kit_replicates_onekit[grepl("_rep", kit_replicates_onekit$kit_id),
                                "wson_conc_field"] - 
            kit_replicates_onekit[!grepl("_rep", kit_replicates_onekit$kit_id), 
                                  "wson_conc_field"]) /
           kit_replicates_onekit[!grepl("_rep", kit_replicates_onekit$kit_id), 
                                 "wson_conc_field"]) * 100)
}

# write out the perc error calcs
#write_csv(npoc_tdn_perc_diff, "./WSOM Extracts Data/wsom_extracts_npoc_tdn_perc_error_calcs.csv")


## 5. Replace Replicate Kit Concentrations w/ Averages -------------------------

# remove replicate rows and cut out unnecessary columns
wsom_extracts_noreps <- wsom_extracts %>%
  filter(!grepl("_rep", kit_id))

# make replacement data frame
wsom_extracts_noreps_replaced <- wsom_extracts_noreps

# replace wsoc and wson concentrations
for (x in 1:nrow(kit_replicates_averages)) {
  
  # set search and replace values
  search_id <- kit_replicates_averages[x,"kit_id"]
  search_location <- kit_replicates_averages[x,"transect_location"]
  new_wsoc_value <- kit_replicates_averages[x,"avg_wsoc_conc"]
  new_wson_value <- kit_replicates_averages[x, "avg_wson_conc"]
  new_wsoc_field_value <- kit_replicates_averages[x, "avg_wsoc_conc_field"]
  new_wson_field_value <- kit_replicates_averages[x, "avg_wson_conc_field"]
  
  for (a in 1:nrow(wsom_extracts_noreps_replaced)) {
    
    # find matches between averaged kits and the larger data frame
    if (wsom_extracts_noreps_replaced[a,"kit_id"] == search_id &
        wsom_extracts_noreps_replaced[a,"transect_location"] 
        == search_location) {
      
      wsom_extracts_noreps_replaced[a,"wsoc_conc"] <- new_wsoc_value
      wsom_extracts_noreps_replaced[a,"wson_conc"] <- new_wson_value
      wsom_extracts_noreps_replaced[a,"wsoc_conc_field"] <- new_wsoc_field_value
      wsom_extracts_noreps_replaced[a,"wson_conc_field"] <- new_wson_field_value
    }
  }
}


## 6. Final Data Edits and Write Out -------------------------------------------

# produce data frame with just ID information and concentrations in
# mg per g dry weight
wsom_conc <- wsom_extracts_noreps_replaced %>%
  select(-c(npoc_mgl,tdn_mgl,water_L,field_soil))

#write_csv(wsom_conc, "./WSOM Extracts Data/wsom_conc_rep_avg.csv")
