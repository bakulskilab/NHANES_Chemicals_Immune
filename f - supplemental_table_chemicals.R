#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
############################################  Supplemental Table 1  ###########################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Purpose: This function organizes the information about the chemicals
#          
# Inputs:   chem_master      - dataframe of all information about the NHANES chemicals
#           conversion       - dataframe of 196 NHANES chemicals in this analysis
#
# Outputs:  supp_table_chems - dataframe of chemicals with all relevant information

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

supplemental_table_chemicals <- function(chem_master,
                                         conversion)
{
  library(tidyverse)
  
  setwd(current_directory)
  
  #TEMPORARY
  # chem_master <- list_master_files$Chemicals
  # conversion <- use_these_chems

  #############################################################################################################
  ########################################## Drop Unnecessary Columns #########################################
  #############################################################################################################
  
  #these are the columns to drop
  keep_cols <- c("chemical_codename",
                 "chemical_codename_use",
                 "chemical_name",
                 "units",
                 "cas_num",
                 "SDDSRVYR",
                 "file_name",
                 "comment_codename",
                 "comments_codename_use",
                 "chem_family",
                 "LOD",
                 "weight_codename")
  
  #drop the columns
  chem_simple <- chem_master %>%
    dplyr::select(all_of(keep_cols))
  
  #remove units from chemical names
  chem_simple$chemical_name <- gsub("\\s\\(([^()]+)\\)$",
                                    "",
                                    chem_simple$chemical_name)

  #############################################################################################################
  ######################################## Keep Only Project Chemicals ########################################
  #############################################################################################################
  
  # Grab the chemicals for this project
  chems <- conversion$chemical_codename_use
  
  # Keep only those chemicals
  chem_subset <- chem_simple %>%
    filter(chemical_codename_use %in% chems)
  
  options(digits=3)
  
  # Keep only some variable columns
  chem_subset_final <- chem_subset %>%
    select(chemical_codename_use,
           chemical_name,
           units,
           cas_num,
           SDDSRVYR,
           file_name,
           chem_family,
           comments_codename_use,
           LOD,
           weight_codename) %>%
    mutate(Matrix = case_when(str_detect(chemical_codename_use, "^LB") ~ "Blood",
                              str_detect(chemical_codename_use, "^UR") ~ "Urine",
                              str_detect(chemical_codename_use, "^SS") ~ "Urine")) %>%
    mutate('Survey Years' = case_when(str_detect(SDDSRVYR, "10") ~ "2017-2018",
                                      str_detect(SDDSRVYR, "9") ~ "2015-2016",
                                      str_detect(SDDSRVYR, "8") ~ "2013-2014",
                                      str_detect(SDDSRVYR, "7") ~ "2011-2012",
                                      str_detect(SDDSRVYR, "6") ~ "2009-2010",
                                      str_detect(SDDSRVYR, "5") ~ "2007-2008",
                                      str_detect(SDDSRVYR, "4") ~ "2005-2006",
                                      str_detect(SDDSRVYR, "3") ~ "2003-2004",
                                      str_detect(SDDSRVYR, "2") ~ "2001-2002",
                                      str_detect(SDDSRVYR, "1") ~ "1999-2000")) %>%
    rename('Chemical Codename (analysis)' = chemical_codename_use,
           'Chemical Name' = chemical_name,
           'Units' = units,
           'CAS number' = cas_num,
           'Survey Cycle' = SDDSRVYR,
           'NHANES File' = file_name,
           'Chemical Family' = chem_family,
           'Comment Codename' = comments_codename_use,
           'LOD' =  LOD,
           'Weight Codename' = weight_codename) %>%
    relocate('Survey Years', .after = 'Survey Cycle')

  
  # Save Table
  print("save chem_master_clean.csv into Table 1 folder")
  setwd(paste0(current_directory, "/Tables - Table 1"))
  write.csv(chem_subset_final, "chem_master_clean.csv", row.names = FALSE)
  
  #############################################################################################################
  #################################### Make a Table of Lipid Adjusted Chems ###################################
  #############################################################################################################
  
  # Subset out the lipid chemicals
  lipids <- conversion %>% filter(str_detect(chemical_name, "(?i)Lipid")) %>% pull(chemical_codename_use)
  lipid_chems <- chem_master %>%
    filter(str_detect(chemical_name, "(?i)Lipid")) %>%
    mutate(chemical_name_only = gsub("\\s\\(([^()]+)\\)$", "", chemical_name)) %>%
    filter(chemical_codename_use %in% lipids) %>%
    select(chemical_codename_use,
           chemical_name_only,
           units,
           cas_num,
           SDDSRVYR,
           file_name,
           chem_family,
           comments_codename_use,
           LOD,
           weight_codename) %>%
    rename('Chemical Codename (analysis)' = chemical_codename_use,
           'Chemical Name' = chemical_name_only,
           'Units' = units,
           'CAS number' = cas_num,
           'Survey Cycle' = SDDSRVYR,
           'NHANES File' = file_name,
           'Chemical Family' = chem_family,
           'Comment Codename' = comments_codename_use,
           'LOD' =  LOD,
           'Weight Codename' = weight_codename)
  
  length(unique(lipid_chems$'Chemical Codename (analysis)'))
  
  #79 lipid chemicals in full list of chemicals
  #56 lipid chemicals in chemicals we're using for this analysis
  
  # Save Table of full lipid chemicals
  print("save chem_master_lipid_clean.csv into Table 1 folder")
  setwd(paste0(current_directory, "/Tables - Table 1"))
  write.csv(lipid_chems, "chem_master_lipid_clean.csv", row.names = FALSE)
  
  
  #############################################################################################################
  #############################################################################################################
  #############################################################################################################
  
  # Reset the directory to the project folder
  setwd(current_directory)
}