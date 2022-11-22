## This script reads in several EXCHANGE data downloads from soil and water
## analyses and pulls them together to make one large datafile for use in
## scripts to look at patterns in WSOC extracts for select EXCHANGE sites and 
## also water-soil relationships across all EXHCANGE sites.
##
## Data are all either downloaded from GDrive or GitHub and pulled into the
## repository associated with this RProject.

## Created: 2022/11/21
## Julia McElhinny

## 1. Load Necessary Packages and Data Files -----------------------------------

require(tidyverse)

# set file location
directory = "./EXCHANGE Downloaded Data"

# pulled from ec1-tai-geochemistry GitHub repo
soils_all <- readRDS(paste0(directory,"/soils_data_merged_withmeta.rds"))

# all from EXCHANGE Data GDrive folder
# metadata for sampling locations and conditions, water NPOC/TDN, 
# TSS, and pH, alkalinity, etc.
sample_metadata <- read.csv(paste0(directory,
                                   "/EC1_Metadata_CollectionLevel.csv"))
water_npoc_tdn <- read.csv(paste0(directory,
                                  "/EC1_Water_NPOC_TDN_L0B_20220601.csv"))
water_tss <- read.csv(paste0(directory,
                             "/EC1_Water_TSS_L0B_20220607.csv"))
water_quality <- read.csv(paste0(directory,
                                 "/EC1_Water_WaterQuality_L0B_20220509.csv"))


## 2. Merge Water Data Together ------------------------------------------------

# make sure transect_location is "Water" not "water"
water_npoc_tdn$transect_location <- str_to_title(water_npoc_tdn$transect_location)

# make list of data frames that need to be merged
water_df_list <- list(water_npoc_tdn, water_tss, water_quality)

# merge data frames
water_all <- water_df_list %>%
  reduce(left_join, by = c("campaign", "kit_id", "transect_location"))

# reduce number of columns
water_all_edit <- water_all[,-c(1,4,7,8,10:13,18:21)]

# edit column names to include "water"
colnames(water_all_edit)[3:9] <- paste("water", colnames(water_all_edit)[3:9],
                                       sep = "_")


## 3. Merge Water Variables with Soil Data Frame -------------------------------

# remove unnecessary columns
soils_all_edit <- soils_all[,-grep("flag", colnames(soils_all))]
soils_all_edit <- soils_all_edit[,-c(7,13,14,17,18,21:25)]

# reorder columns to bring some metadata next to kit_id
soils_all_edit <- soils_all_edit[,c(1,2,14,15,3:13,16:18)]

water_soil_characteristics <- inner_join(water_all_edit, soils_all_edit, 
                                         by = c("campaign", "kit_id"))

## 4. Add Collection Metadata --------------------------------------------------

# pull lat long info to make easy data frame for all GIS work
collection_coordinates <- sample_metadata %>%
  select(contains(c("lat", "long")))

collection_coordinates$kit_id <- sample_metadata$kit_id
collection_coordinates <- collection_coordinates[,c(11,1:10)]

# pivot the columns into rows for latitude
collection_coordinates <- collection_coordinates %>%
  pivot_longer(c(2:6), names_to = "transect_location", values_to = "latitude")

# mutate transect_location
locations <- as.data.frame(unique(collection_coordinates$transect_location))
colnames(locations) <- "old"

locations$new <- c("Water", "Sediment", "Wetland", "Transition", "Upland")

for (n in 1:nrow(locations)) {
  search <- locations[n,1]
  replace_value <- locations[n,2]
  
  # find matching values to "search" and replace with "replace"
  collection_coordinates$transect_location <- sapply(
                                      collection_coordinates$transect_location,
                                                     function(x){x[x==search] 
                                                       <- replace_value
                                                     
                                                     return(x)
                                                     })
}

# change column names of long values to match transect location
colnames(collection_coordinates)[2:6] <- locations$new

# pivot remaining columns into rows for longitude
collection_coordinates <- collection_coordinates %>%
  pivot_longer(c(2:6), names_to = "transect_location", values_to = "longitude")


# make metadata into tidyish format to merge with other parameters
# take water info out

