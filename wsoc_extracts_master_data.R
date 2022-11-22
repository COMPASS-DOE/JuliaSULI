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
require(data.table) # for editing to tidy data
require(stringi) # for editing text strings

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


## 4. Add Collection Metadata to Master Data Set -------------------------------

# pull lat long info to make easy data frame for all GIS work
collection_coordinates <- sample_metadata %>%
  select(contains(c("lat", "long")))

collection_coordinates$kit_id <- sample_metadata$kit_id
collection_coordinates <- collection_coordinates[,c(11,1:10)]

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

# make sure kit_id and transect_location are factor variables
collection_coordinates[,c("kit_id", "transect_location")] <- lapply(
  collection_coordinates[,c("kit_id", "transect_location")],
  function(x) as.factor(x)
)

# write out rds file to use in GIS-related scripts
write_rds(collection_coordinates, 
          "./Processed Data/sample_collection_coordinates.rds")

# add water info
water_metadata <- sample_metadata[,c(1,6,7)]

water_soil_characteristics <- left_join(water_soil_characteristics, 
                                        water_metadata,
                                        by = "kit_id")


## 5. Add Final Edits to the Data Frame ----------------------------------------

# add water body variable
water_soil_characteristics <- water_soil_characteristics %>%
  mutate(water_body = case_when(
    region == "Great Lakes" ~ "Great Lakes",
    state == "VA" | state == "MD" ~ "Chesapeake Bay",
    state == "PA" | state == "DE" | state == "NJ" ~ "Delaware Bay"
  ))

# edit water_systemtype variable
water_soil_characteristics$water_systemtype <- stri_replace(
  water_soil_characteristics$water_systemtype,
  fixed = ", ", "/"
)

# make variables factor variables
water_soil_characteristics[,c("kit_id", "transect_location", "state",
                              "region", "water_body", "water_systemtype",
                              "soil_horizon", "visible_iron_oxidation",
                              "visible_minerals", "visible_white_flakes")] <-
  lapply(water_soil_characteristics[,c("kit_id", "transect_location", "state",
                                       "region", "water_body", 
                                       "water_systemtype", "soil_horizon", 
                                       "visible_iron_oxidation", 
                                       "visible_minerals", 
                                       "visible_white_flakes")],
         function(x) as.factor(x)
  )

# reorder variables
water_soil_characteristics <- water_soil_characteristics[,c(1,2,10,11,28,27,3:9,
                                                            26,12:25)]

# write out data frame as rds file
write_rds(water_soil_characteristics, 
          "./Processed Data/EC1_water_soil_data.rds")
