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

require(pacman)
pacman::p_load(tidyverse,
               data.table, # for editing to tidy data
               stringi) # for editing text strings

# set file location
directory = "./Downloaded Data"

# pulled from #ec1-tai-geochemistry GitHub repo
soils_all <- readRDS(paste0(directory,"/soils_data_merged_withmeta.rds"))

# all from EXCHANGE Data GDrive folder
# metadata for sampling locations and conditions, water NPOC/TDN, 
# TSS, and pH, alkalinity, etc.
sample_metadata <- read.csv(paste0(directory,
                                   "/EC1_Metadata_CollectionLevel.csv"))
kit_metadata <- read.csv(paste0(directory,
                                "/EC1_Metadata_KitLevel.csv"))
water_npoc_tdn <- read.csv(paste0(directory,
                                  "/EC1_Water_NPOC_TDN_L0B_20220601.csv"))
water_tss <- read.csv(paste0(directory,
                             "/EC1_Water_TSS_L0B_20220607.csv"))
water_quality <- read.csv(paste0(directory,
                                 "/EC1_Water_WaterQuality_L0B_20220509.csv"))
water_spectral_indices <- read.csv(paste0(directory, 
                                          "/EC1_Water_CDOM_SpectralIndices_20221011.csv"))
water_peak_picks <- read.csv(paste0(directory,
                                    "/EC1_Water_CDOM_DOCnormalizedpeakpicks_20221011.csv"))



## 2. Merge Water Data Together ------------------------------------------------

# make sure transect_location is "Water" not "water"
water_npoc_tdn$transect_location <- str_to_title(water_npoc_tdn$transect_location)

# change column name with kit info to "kit_id" for CDOM data
water_spectral_indices <- water_spectral_indices %>%
  rename(kit_id = Sample_ID) %>% 
  mutate(campaign = "EC1",
         transect_location = "Water")
water_peak_picks <- water_peak_picks %>%
  rename(kit_id = Sample_ID) %>%
  mutate(campaign = "EC1",
         transect_location = "Water")

# make list of data frames that need to be merged
water_df_list <- list(water_npoc_tdn, water_tss, water_quality,
                      water_spectral_indices, water_peak_picks)

# merge data frames
water_all <- water_df_list %>%
  reduce(left_join, by = c("campaign", "kit_id", "transect_location")) %>%
  select(-Sample_Description, -date)

# edit column names to include "water"
colnames(water_all)[3:44] <- paste("water", colnames(water_all)[3:44],
                                       sep = "_")

# merge with kit-level metadata
water_all <- full_join(water_all, kit_metadata,
                       by = "kit_id")

# merge with collection-level metadata that has to do with water
water_metadata <- sample_metadata %>%
  select(kit_id, contains("water"))

water_all <- full_join(water_all, water_metadata) %>%
  mutate(campaign = "EC1")

# write out RDS for all water data
#write_rds(water_all, "./water_data_mergedwithallmeta.rds")


## 3. Merge Water Variables with Soil Data Frame -------------------------------

# reorder columns to bring some metadata next to kit_id
soils_all_edit <- soils_all[,c(1:3,24:28,4:23,29:35)]

# work w/ K046 - transition measurements are supposed to be wetland (mess-up in labelling)
kit46 <- soils_all_edit %>%
  filter(kit_id == "K046") %>%
  filter(transect_location == "Wetland"|transect_location == "Transition")

for (a in 1:ncol(kit46)) {
  # if a cell is NA for the wetland row, replace it with the value from the transition row
   if (is.na(kit46[1,a])) {
     kit46[1,a] <- kit46[2,a]
   }
}

kit46 <- kit46 %>% filter(transect_location == "Wetland")

soils_all_edit <- soils_all_edit %>%
  filter(!(kit_id == "K046" & transect_location == "Wetland")) %>%
  filter(!(kit_id == "K046" & transect_location == "Transtion")) %>%
  rbind(kit46)
  
# play with collection level metadata for all the soils
soils_metadata <- sample_metadata %>%
  select(kit_id, contains("sediment")|contains("wetland")|contains("transition")|contains("upland")) %>%
  pivot_longer(-kit_id,
               names_to = c("transect_location", ".value"),
               names_pattern = "^([A-Za-z]+)_([A-Za-z]+_?[A-Za-z]+_?[A-Za-z]+)$") %>%
  mutate(transect_location = str_to_title(transect_location))
# this names_pattern line basically separates the character chunks at the first underscore ONLY
# [A-Za-z] is looking to match alphabetic characters
# all the ? after the first underscore means those characteristics are optional (or don't show up in all the column names)
# other names_pattern tried: "(.+)_(.+)" split at every underscore, "(.+)_(.+_?.+_?.+)" worked on some of the column names but not all

# merge this metadata with the soils info
soils_all_edit <- full_join(soils_all_edit, soils_metadata,
                       by = c("kit_id", "transect_location"))

# some issues were created with the merging - added kits to the list that weren't collected
# make a list of the kits missing lat and long information
soils_not_collected <- soils_all_edit %>%
  filter(is.na(latitude)) %>%
  filter(kit_id != "K014") %>% #K014 is missing metadata - was collected!!
  mutate(sample_ID = paste(kit_id, transect_location, sep = "_"))

soils_all_final <- soils_all_edit %>%
  select(-soil_horizon.x) %>%
  rename(soil_horizon = soil_horizon.y) %>% # the second soil horizon has more rows identified 
  mutate(campaign = "EC1") %>% # make sure this column is filled in for all rows
  mutate(sample_ID = paste(kit_id, transect_location, sep = "_")) %>%
  filter(!(sample_ID %in% soils_not_collected$sample_ID)) %>%
  select(-sample_ID)

# write out soil .rds file with updated metadata merge
#write_rds(soils_all_final, "./soils_data_mergedwithallmeta.rds")
  
# make data frame with water and soil data together and write out as .rds file
water_soil <- full_join(water_all, soils_all_final)


## 4. Add Collection Metadata to Master Data Set -------------------------------

# pull lat long info to make easy data frame for all GIS work
collection_coordinates <- sample_metadata %>%
  select(kit_id, contains(c("lat", "long")))

# pivot the columns into rows for latitude (using data.table)
# pivot_longer function in tidyr package was difficult to break down two 
# different groups of columns into one "key" column
dt <- data.table(collection_coordinates)

dt_edit <- melt(dt, id.vars = "kit_id",
                measure.vars = patterns(latitude = "_latitude",
                                        longitude = "_longitude"),
                variable.name = "transect_location")

# make data table back into data frame
collection_coordinates <- as.data.frame(dt_edit)

# change names of transect locations
collection_coordinates <- collection_coordinates %>%
  mutate(transect_location = case_when(
              transect_location == 1 ~ "Water",
              transect_location == 2 ~ "Sediment",
              transect_location == 3 ~ "Wetland",
              transect_location == 4 ~ "Transition",
              transect_location == 5 ~ "Upland"))

# another way to do the above steps using pivot_longer instead of melt():
#collection_coordinates <- collection_coordinates %>%
 # dplyr::select(kit_id, contains("latitude")|contains("longitude")) %>%
  #pivot_longer(-kit_id, 
   #            names_to = c("transect_location", ".value"),
    #           names_pattern = "(.+)_(.+)") %>%
  #mutate(transect_location = str_to_title(transect_location))

# make sure kit_id and transect_location are factor variables
collection_coordinates[,c("kit_id", "transect_location")] <- lapply(
  collection_coordinates[,c("kit_id", "transect_location")],
  function(x) as.factor(x)
)

# write out a csv file to use in other scripts
#write_csv(collection_coordinates, "./kit_sample_locations_lat_long.csv")


## 5. Add Final Edits to the Data Frame ----------------------------------------

# add water body variable
water_soil <- water_soil %>%
  mutate(water_body = case_when(
    region == "Great Lakes" ~ "Great Lakes",
    state == "VA" | state == "MD" ~ "Chesapeake Bay",
    state == "PA" | state == "DE" | state == "NJ" ~ "Delaware Bay"
  )) %>%
  relocate(water_body, .after = region)

# edit water_systemtype variable
water_soil$water_systemtype <- stri_replace(
  water_soil$water_systemtype,
  fixed = ", ", "/"
)

# make variables factor variables
water_soil[,c("kit_id", "transect_location", "state",
                              "region", "water_body", "water_systemtype",
                              "soil_horizon", "visible_iron_oxidation",
                              "visible_minerals", "visible_white_flakes")] <-
  lapply(water_soil[,c("kit_id", "transect_location", "state",
                                       "region", "water_body", 
                                       "water_systemtype", "soil_horizon", 
                                       "visible_iron_oxidation", 
                                       "visible_minerals", 
                                       "visible_white_flakes")],
         function(x) as.factor(x)
  )

# reorder variables
water_soil <- water_soil %>%
  relocate(c(collection_date:samples_collected), .after = kit_id)
  

# write out data frame as rds file
#write_rds(water_soil, "./EC1_water_soil_data.rds")
