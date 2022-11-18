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

### QUESTION for team: We are setting blanks to 0 when the blank is less than LOD
### Does that make sense? Setting blanks to 0 when they are not? 


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
#lod_npoc <- 0.076
#lod_tdn <- 0.014

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

# make the date column numerical
npoc_raw$date <- as.numeric(npoc_raw$date)

# mark rows that are duplicates
npoc_raw$dups <- duplicated(npoc_raw$sample_name)

## replace values over cal curve with the rerun NPOC values
# filter out the duplicate rows into another data frame
npoc_raw_dups <- npoc_raw %>%
  filter(dups == TRUE)

# make a dataset without the duplicates
npoc_raw_noduplicates <- npoc_raw %>%
  filter(dups == FALSE)

# make a dataset to show the replacements
npoc_raw_replaced <- npoc_raw_noduplicates

# loop through the rows needing replacement, check conditions, and replace values
for (n in 1:nrow(npoc_raw_dups)) {
  
  # find the rows in the original dataset that match names needing to be replaced
  for (i in 1:nrow(npoc_raw_replaced)) {
   if (npoc_raw_replaced$sample_name[i] == npoc_raw_dups$sample_name[n]) {
     
     # check if the date of the duplicates is greater than the date of the original data 
     # assumption here that rerun data will be more recent and will be the desired values to process (because they were rerun for a reason)
     if (npoc_raw_dups$date[n] > npoc_raw_replaced$date[i]) {
       
       # replace values
       npoc_raw_replaced$npoc_raw[i] <- npoc_raw_dups$npoc_raw[n]
       
       # replace the date of measurement
       npoc_raw_replaced$date[i] <- npoc_raw_dups$date[n]
     }
   }
  }
}

# 4. Calculate blanks and add to data ------------------------------------------

blanks <- npoc_raw %>% 
  filter(grepl("^Blank", sample_name)) %>% 
  group_by(date) %>% 
  summarize(npoc_blank_raw = round(mean(npoc_raw[!is.na(npoc_raw)]), 2), 
            tdn_blank_raw = round(mean(tdn_raw[!is.na(tdn_raw)]), 2)) %>% 
  mutate(npoc_blank = ifelse(npoc_blank_raw > lod_npoc, npoc_blank_raw, 0), 
         tdn_blank = ifelse(tdn_blank_raw > lod_tdn, tdn_blank_raw, 0)) %>% 
  select(date, npoc_blank, tdn_blank)


# 5. Add blanks data -----------------------------------------------------------

npoc_blank_corrected <- npoc_raw %>% 
  filter(grepl("EC1_K", sample_name)) %>% # filter to EC1 samples only
  mutate(campaign = "EC1", 
         kit_id = substr(sample_name, 5, 9), 
         transect_location = "water") %>% 
  inner_join(blanks, by = "date") %>% 
  mutate(npoc_mgl = npoc_raw - npoc_blank, 
         tdn_mgl = tdn_raw - tdn_blank)


# 6. Clean data ----------------------------------------------------------------

## Helper function to calculate mean if numeric, otherwise first (needed to 
## preserve dates, which are different for duplicated kits)
mean_if_numeric <- function(x){
  ifelse(is.numeric(x), mean(x, na.rm = TRUE), first(x))
}

## Another step before finalizing is taking care of pesky duplicates from reruns
npoc_duplicates_removed <- npoc_blank_corrected %>% 
  select(campaign, transect_location, kit_id, date, npoc_mgl, tdn_mgl, npoc_blank, tdn_blank) %>% 
  group_by(kit_id) %>% 
  summarize(across(everything(), .f = mean_if_numeric))

## The last step is flagging data
npoc_raw_flags <- npoc_duplicates_removed %>% 
  ## First, round each parameter to proper significant figures
  mutate(npoc_mgl = round(npoc_mgl, 2), 
         tdn_mgl = round(tdn_mgl, 3)) %>% 
  ## Second, add flags for outside LOD
  mutate(npoc_flag = ifelse(npoc_mgl < lod_npoc | npoc_mgl > 30, "npoc outside range", NA), #per cal curve upper limit
         tdn_flag = ifelse(tdn_mgl < lod_tdn | tdn_mgl > 3, "tdn outside range", NA))

npoc <- npoc_raw_flags %>% 
  select(date, campaign, kit_id, transect_location, npoc_mgl, tdn_mgl, contains("_flag"))


# 7. Write data ----------------------------------------------------------------
date_updated <- "20220601"

write_csv(npoc, paste0("Data/Processed/EC1_Water_NPOC_TDN_L0B_", date_updated, ".csv"))

