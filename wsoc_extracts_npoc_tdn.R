## This script imports raw data for NPOC and TDN measured using a Shimadzu TOC-L
## at PNNL MCRL and exports clean, Level 0B QC'ed data. 
## Data are read in from GitHub
## 
## Created: 2022-01-15 (Updated 2022-03-30 with OO)
## Peter Regier
## Edited: 2022-11 by Julia McElhinny for use with WSOC Extracts Data
## w/ code chunks from AMP edited version used for TEMPEST porewater data
##
# #############
# #############

# 1. Setup ---------------------------------------------------------------------

# load packages
require(pacman)
pacman::p_load(tidyverse, # keep things tidy
               janitor, # useful for simplifying column names
               googlesheets4, # read_sheet 
               googledrive) # drive_ functions

## Set theme
theme_set(theme_bw())

## Set LOD (for all data after 10/11/2021)
## If any data were before 10/11/21, lod_npoc = 0.27, lod_tdn = 0.070
lod_npoc <- 0.076
lod_tdn <- 0.014

## Set GitHub repository folder for raw data files
directory = "./Raw NPOC Data"

# 2. Functions -----------------------------------------------------------------

## Create a function to read in data
read_data <- function(data){
  # First, scrape date from filename
  # pull characters from filename not filepath/directory name
  date <- substr(data, nchar(directory) + 2, nchar(directory) + 9) 
  # Second, read in data
  read_delim(file = data, skip = 10, delim = "\t") %>% 
    rename(sample_name = `Sample Name`, 
           npoc_raw = `Result(NPOC)`, 
           tdn_raw = `Result(TN)`) %>% 
    select(sample_name, npoc_raw, tdn_raw) %>% 
    mutate(date = date)
}

## Build function to read in ReadMes
read_mes <- function(readme){
  # First, scrape date from filename
  date <- substr(readme, nchar(directory) + 2, nchar(directory) + 9)
  # Second, read in Read Me
  readxl::read_excel(path = readme, sheet = 1) %>% 
    rename(sample_name = `Sample Name`) %>% 
    select(sample_name, Action, Dilution) %>% 
    mutate(date = date)
}

# 3. Import data ---------------------------------------------------------------

## Create a list of files to download
files <- list.files(path = directory, pattern = "Summary", full.names = TRUE) 
ReadMes <- list.files(path = directory, pattern = "Readme", full.names = TRUE) 

## Read in NPOC data, filter to EC1 samples
npoc_raw <- files %>% 
  map_df(read_data) %>% 
  filter(grepl("EC1", sample_name)) %>% 
  bind_rows()

## Read in ReadMes, filter to EC1 samples only
readmes_all <- ReadMes %>% 
  map_df(read_mes) %>% 
  filter(grepl("EC1", sample_name)) %>%
  bind_rows() 

# 4. Work with Duplicate NPOC Values -------------------------------------------

## Replace values over cal curve with the rerun NPOC values and the dilution 
## factors
# make the date column numerical
npoc_raw$date <- as.numeric(npoc_raw$date)
readmes_all$date <- as.numeric(readmes_all$date)

# mark rows that are duplicates
npoc_raw$dups <- duplicated(npoc_raw$sample_name)
readmes_all$dups <- duplicated(npoc_raw$sample_name)

# filter out the duplicate rows into another data frame
npoc_raw_dups <- npoc_raw %>%
  filter(dups == TRUE)
readmes_all_dups <- readmes_all %>%
  filter(dups == TRUE)

# combine the duplicate data frames
npoc_raw_dups <- left_join(npoc_raw_dups, readmes_all_dups, 
                           by = c("sample_name", "date", "dups"))

# make a dataset without the duplicates
npoc_raw_noduplicates <- npoc_raw %>%
  filter(dups == FALSE)
readmes_all_noduplicates <- readmes_all %>%
  filter(dups == FALSE)

# make a dataset to show the replacements
npoc_raw_replaced <- npoc_raw_noduplicates
readmes_all_dilution_replaced <- readmes_all_noduplicates

# edit date column to have a date for NPOC measurement and one for TDN
npoc_raw_replaced$npoc_date <- npoc_raw_replaced$date
colnames(npoc_raw_replaced)[4] <- "tdn_date"

# loop through the rows needing replacement, check conditions, replace values
for (n in 1:nrow(npoc_raw_dups)) {
  
  # find the rows in the original dataset that match names that need replaced
  for (i in 1:nrow(npoc_raw_replaced)) {
   if (npoc_raw_replaced$sample_name[i] == npoc_raw_dups$sample_name[n]) {
     
     # check if the date of the duplicates is greater than the date of the 
     # original data 
     # assumption here that rerun data will be more recent and will be the 
     # desired values to process (because they were rerun for a reason)
     if (npoc_raw_dups$date[n] > npoc_raw_replaced$npoc_date[i]) {
       
       # replace values
       npoc_raw_replaced$npoc_raw[i] <- npoc_raw_dups$npoc_raw[n]
       
       # replace the date of measurement
       npoc_raw_replaced$npoc_date[i] <- npoc_raw_dups$date[n]
     }
   }
  }
}

# edit dilution column to have one factor for TDN and one for NPOC
readmes_all_dilution_replaced$npoc_dilution <- 
  readmes_all_dilution_replaced$Dilution
colnames(readmes_all_dilution_replaced)[3] <- "tdn_dilution"

# same loop, but for dilution factors
for (n in 1:nrow(npoc_raw_dups)) {
  
  # find the rows in the original dataset that match names that need replaced
  for (i in 1:nrow(readmes_all_dilution_replaced)) {
    if (readmes_all_dilution_replaced$sample_name[i] == 
        npoc_raw_dups$sample_name[n]) {
      
      # check if the date of the duplicates is greater than the date of the 
      # original data 
      # assumption here that rerun data will be more recent and will be the 
      # desired values to process (because they were rerun for a reason)
      if (npoc_raw_dups$date[n] > readmes_all_dilution_replaced$date[i]) {
        
        # replace values
        readmes_all_dilution_replaced$npoc_dilution[i] <- 
          npoc_raw_dups$Dilution[n]
        
        # replace the date of measurement
        readmes_all_dilution_replaced$date[i] <- npoc_raw_dups$date[n]
      }
    }
  }
}

# 5. Calculate process blanks --------------------------------------------------

# pull process blanks out of dataset
# calculate NPOC average and TDN if undiluted
# make sure blanks are either 0 or higher than instrument LOD
blanks <- npoc_raw_replaced %>% 
  filter(grepl("blank-filter", sample_name)) %>% 
  group_by(npoc_date, tdn_date) %>%
  summarize(npoc_blank_raw = round(mean(npoc_raw[!is.na(npoc_raw)]), 2),
            tdn_blank_raw = round(mean(tdn_raw[!is.na(tdn_raw)]), 2)) %>% 
  mutate(npoc_blank = ifelse(npoc_blank_raw > lod_npoc, npoc_blank_raw, 0),
         tdn_blank = ifelse(tdn_blank_raw > lod_tdn, tdn_blank_raw, 0)) %>% 
  select(npoc_date, npoc_blank, tdn_date, tdn_blank)

colnames(blanks)[1] <- "npoc_blank_date"

# 6. Blank Correction ----------------------------------------------------------

npoc_blank_corrected <- npoc_raw_replaced %>% 
  filter(grepl("EC1_K", sample_name)) %>% # filter to EC1 samples only
  mutate(campaign = "EC1", 
         kit_id = substr(sample_name, 5, 8), 
         transect_location = case_when(
           substr(sample_name, 10, 11) == "UP" ~ "Upland",
           substr(sample_name, 10, 10) == "T" ~ "Transition",
           substr(sample_name, 10, 10)  == "W" ~ "Wetland")) %>% 
  inner_join(blanks, by = "tdn_date") %>% 
  mutate(npoc_bc = npoc_raw - npoc_blank, 
         tdn_bc = tdn_raw - tdn_blank)

# add blank-sol and blank-filter values back in for CDOM normalizations
npoc_blank_corrected <- bind_rows(npoc_blank_corrected, 
                              filter(npoc_raw_replaced, 
                                     grepl("EC1_b", sample_name)))


# 7. Dilution Correction -------------------------------------------------------

colnames(readmes_all_dilution_replaced)[4] <- "npoc_date"

samples_dilution_corrected <- npoc_blank_corrected %>%
  left_join(readmes_all_dilution_replaced, by = c("sample_name", "npoc_date", 
                                                  "dups")) %>% 
  mutate(npoc_mgl = ifelse(!is.na(npoc_bc), npoc_bc * npoc_dilution,
                           npoc_raw * npoc_dilution),
         tdn_mgl = ifelse(!is.na(tdn_bc), tdn_bc * tdn_dilution, 
                          tdn_raw * tdn_dilution), 
         npoc_mgl = as.numeric(npoc_mgl), npoc_mgl = round(npoc_mgl, 2),
         tdn_mgl= as.numeric(tdn_mgl), tdn_mgl= round(tdn_mgl, 2))

# 8. Final Dataset Prep --------------------------------------------------------

# overwrite blank-sol TN values because they were diluted MQ which doesn't 
# make sense; values won't be used in analysis anyways
samples_dilution_corrected[c(35,37,39), "tdn_mgl"] <- 9999999

# select desired columns
wsoc <- samples_dilution_corrected %>%
  select(sample_name, tdn_date, npoc_date, campaign, kit_id, transect_location, 
         npoc_mgl, tdn_mgl)

# 9. Write out data ------------------------------------------------------------
date_updated <- "20221118"

write_csv(wsoc, paste0("./Processed Data/EC1_WSOC_Extracts_NPOC_TDN_L0B_", 
                       date_updated, ".csv"))
