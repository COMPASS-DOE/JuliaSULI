## This script imports raw data for NPOC and TDN measured using a Shimadzu TOC-L
## at PNNL MCRL and exports clean, Level 0B QC'ed data. 
## Data are read in from GitHub
## 
## Created: 2022-01-15 (Updated 2022-03-30 with OO)
## Peter Regier
## Edited: 2022-11 by Julia McElhinny for use with WSOM Extracts Data
## w/ code chunks from AMP's TEMPEST porewater data processing script
##
# #############
# #############

# 1. Setup ---------------------------------------------------------------------

# load packages
require(pacman)
pacman::p_load(tidyverse, # keep things tidy
               janitor, # useful for simplifying column names
               readxl) # read in excel spreadsheets

## Set theme
theme_set(theme_bw())

## Set GitHub repository folder for raw data files
directory = "./WSOM Extracts Data/raw npoc"

# 2. Functions -----------------------------------------------------------------

## Create a function to read in data
read_data <- function(data){
  # First, scrape date from filename
  # pull characters from filename not filepath/directory name
  # when the files are in the GitHub repo, the filename includes the entire path...
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

## Read in NPOC data
npoc_raw <- files %>% 
  map_df(read_data) %>%
  bind_rows()

## Read in ReadMes, filter to EC1 samples only
readmes_all <- ReadMes %>% 
  map_df(read_mes) %>% 
  filter(grepl("EC1_K", sample_name)) %>%
  bind_rows() 

## Read in LOD information for TOC instrument
LOD <- read_excel("TOC_MCRL_LOD.xlsx")


# 4. Work with LOD -------------------------------------------------------------

## edit column names of LOD data frame for easier reference
colnames(LOD) <- c("Date_LOD_run", "LOD_Start_Date", "LOD_End_Date", "LOD_NPOC", "LOD_TN", 
                   "LOD_DIC", "LOD_TSS", "Instrument_Notes")

# fill NA end date with System date
LOD[is.na(LOD$LOD_End_Date),3] <- as.numeric(gsub("-", "", Sys.Date()))

## make date columns numeric for math later
LOD$Date_LOD_run <- as.numeric(LOD$Date_LOD_run)
npoc_raw$date <- as.numeric(npoc_raw$date)

## create a blank list to populate
raw_data_LOD_list <- list()

## for loop to determine which LOD the instrument runs fall within
## this is done by calculating which LOD period start date is closest to the date of the instrument run
for (a in 1:length(unique(npoc_raw$date))) {
  
  ## pull focus date
  focus_date <- unique(npoc_raw$date)[a]
  
  ## filter LOD dataset to just the dates that are less than the desired date
  ## also calculate difference between each of the LOD start dates and the desired date
  target_LOD <- LOD %>%
    filter(Date_LOD_run < focus_date) %>%
    mutate(difference = Date_LOD_run - focus_date) %>%
    filter(difference == max(difference))
  
  ## pull the LOD run start date that is closest to the instrument run date
  target_LOD_date <- target_LOD$Date_LOD_run[1]
  
  ## filter the blanks dataset to all the ones with the focus date
  ## and then add the LOD date as a new column for later merging
  npoc_raw_LOD <- npoc_raw %>%
    filter(date == focus_date) %>%
    mutate(LOD_date = target_LOD_date)
  
  ## populate the empty list
  raw_data_LOD_list[[a]] <- npoc_raw_LOD
  
}

## bind all rows together
raw_data_LOD <- bind_rows(raw_data_LOD_list)

## merge TN and NPOC LOD values to edit blanks that fall below these cutoffs
raw_data_LOD <- left_join(raw_data_LOD, select(LOD, Date_LOD_run, LOD_NPOC, LOD_TN),
                          by = c("LOD_date" = "Date_LOD_run"))


# 5. More Data Set-Up ----------------------------------------------------------

## Filter raw data to just EXCHANGE samples based on "EC1_K" name pattern
samples <- raw_data_LOD %>% 
  filter(grepl("EC1_K", sample_name))

## create a data frame with process blanks
process_blanks <- raw_data_LOD %>%
  filter(grepl("EC1_blank", sample_name))

## Pull information about the calibration curve in another data frame
cal_curve <- raw_data_LOD %>%
  filter(grepl("STD_|STD-", sample_name)) %>%
  group_by(date) %>%
  # pull the value of the top point of the calibration curve for TDN and for NPOC by instrument run date
  # pulling the max value because the cal curve could have been rerun at the end to handle values over it
  summarize(npoc_calib_upper_limit = max(npoc_raw, na.rm = TRUE),
            tdn_calib_upper_limit = max(tdn_raw, na.rm = TRUE)) %>%
  # make sure an NA is in the data frame for cal curves that were not run on specific day
  mutate(npoc_calib_upper_limit = ifelse(npoc_calib_upper_limit > 0,
                                         npoc_calib_upper_limit,
                                         NA),
         tdn_calib_upper_limit = ifelse(tdn_calib_upper_limit > 0,
                                        tdn_calib_upper_limit,
                                        NA)) 


# 6. Calculate Process Blanks --------------------------------------------------

## join readmes with process blanks to determine which process blanks were run undiluted
blanks_dilution <- process_blanks %>%
  mutate(date = as.character(date)) %>%
  left_join(select(readmes_all, sample_name, Dilution, date), by = c("sample_name", "date")) %>%
  # mark dilution of 2 for NPOC as NA - don't want to use these values
  mutate(npoc_raw = ifelse(Dilution > 1,
                           NA,
                           npoc_raw))

## calculate blanks via average
blanks_final <- blanks_dilution %>%
  # filter to just blank-filter blanks
  filter(grepl("filter", sample_name)) %>%
  # check against LOD - should not be an issue
  mutate(npoc_raw = ifelse(!is.na(npoc_raw),
                           ifelse(npoc_raw < LOD_NPOC,
                                  LOD_NPOC,
                                  npoc_raw),
                           NA),
         tdn_raw = ifelse(!is.na(tdn_raw),
                          ifelse(tdn_raw < LOD_TN,
                                 LOD_TN,
                                 tdn_raw),
                          NA)) %>%
  # calculate averages
  summarize(npoc_blank = mean(npoc_raw, na.rm = TRUE),
            tdn_blank = mean(tdn_raw, na.rm = TRUE))

## IMPORTANT NOTE - TDN blank may need to be rethought
## process blanks run 2x diluted only for TDN but wasn't completely a 2x dilution
## ... because process blanks are already mostly MQ water


# 7. Flag Duplicate Samples (and merge cal curve info) -------------------------

## for these sample runs, samples were rerun rediluted if they were outside cal curve
## flag these out of cal curve samples to be disregarded later
samples_cal_curve <- samples %>%
  left_join(cal_curve, by = "date") %>%
  mutate(npoc_cal_flag = ifelse(npoc_raw > npoc_calib_upper_limit,
                                "above cal curve",
                                NA),
         tdn_cal_flag = ifelse(tdn_raw > tdn_calib_upper_limit,
                               "above cal curve",
                               NA)) %>%
  # also flag any samples below LOD - shouldn't be an issue but want to check
  mutate(npoc_LOD_flag = ifelse(npoc_raw < LOD_NPOC,
                                "below LOD",
                                NA),
         tdn_LOD_flag = ifelse(tdn_raw < LOD_TN,
                               "below LOD",
                               NA)) %>%
  # combine flags into one flag column
  mutate(npoc_flag = ifelse(!is.na(npoc_cal_flag),
                            npoc_cal_flag,
                            ifelse(!is.na(npoc_LOD_flag),
                                   npoc_LOD_flag,
                                   NA)),
         tdn_flag = ifelse(!is.na(tdn_cal_flag),
                            tdn_cal_flag,
                            ifelse(!is.na(tdn_LOD_flag),
                                   tdn_LOD_flag,
                                   NA))) %>%
  select(-npoc_cal_flag, -tdn_cal_flag, -npoc_LOD_flag, -tdn_LOD_flag)


# 8. Blank Correct and Dilution Correct Samples --------------------------------

## subtract blank values from npoc and tdn values
samples_bc <- samples_cal_curve %>%
  mutate(npoc_bc = npoc_raw - blanks_final$npoc_blank,
         tdn_bc = tdn_raw - blanks_final$tdn_blank)

## merge readme info in to dilution correct values
samples_bc_dc <- samples_bc %>%
  mutate(date = as.character(date)) %>%
  left_join(readmes_all, by = c("sample_name", "date")) %>%
  mutate(npoc_mgl = npoc_bc * Dilution,
         tdn_mgl = tdn_bc * Dilution, 
         npoc_mgl = as.numeric(npoc_mgl), 
         npoc_mgl = round(npoc_mgl, 2),
         tdn_mgl= as.numeric(tdn_mgl), 
         tdn_mgl= round(tdn_mgl, 2))

## replace values above the cal curve with NA to indicate replacement is necessary
samples_bc_dc <- samples_bc_dc %>%
  mutate(npoc_mgl_edit = ifelse(!is.na(npoc_flag),
                           NA,
                           paste0(npoc_mgl))) %>%
  # cut down columns for final edits
  select(sample_name, date, npoc_mgl, npoc_mgl_edit, tdn_mgl, npoc_flag, tdn_flag)


# 9. Final Data Edits and Write Out --------------------------------------------

## edit blanks for binding back to other samples
blanks <- blanks_dilution %>%
  group_by(sample_name) %>%
  summarize(npoc_mgl = sum(npoc_raw, na.rm = TRUE),
            tdn_mgl = sum(tdn_raw, na.rm = TRUE))

## deal with duplicates (replace better npoc values), extract kit_id and transect_location
wsom <- samples_bc_dc %>%
  group_by(sample_name) %>%
  mutate(npoc_mgl_edit = as.numeric(npoc_mgl_edit)) %>%
  summarize(npoc_mgl = sum(npoc_mgl_edit, na.rm = TRUE),
            tdn_mgl = sum(tdn_mgl, na.rm = TRUE)) %>%
  mutate(kit_id = substr(sample_name, 5, 8), 
         transect_location = case_when(
           substr(sample_name, 10, 11) == "UP" ~ "Upland",
           substr(sample_name, 10, 10) == "T" ~ "Transition",
           substr(sample_name, 10, 10)  == "W" ~ "Wetland")) %>%
  relocate(c(kit_id, transect_location), .before = sample_name)

## add blanks to final data frame for writing out
wsom_final <- bind_rows(wsom, blanks) %>%
  mutate(campaign = "EC1") %>%
  relocate(campaign, .before = kit_id) %>%
  # replace TDN value for blank solution since they were run diluted - doesn't quite make sense
  mutate(tdn_mgl = ifelse(grepl("blank-sol", sample_name),
                          NA,
                          tdn_mgl))

## final write out
#write_csv(wsom_final, "./WSOM Extracts Data/EC1_WSOM_Extracts_NPOC_TDN_L0B.csv")
  
