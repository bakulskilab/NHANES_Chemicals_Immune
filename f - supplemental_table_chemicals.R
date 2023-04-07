#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
############################################  Supplemental Table 1  ###########################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Purpose: This function organizes the information about the chemicals
#          
# Inputs:   chem_master      - dataframe of all information about the NHANES chemicals
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
  
  # Save Table
  print("save chem_master_clean.csv into Table 1 folder")
  setwd(paste0(current_directory, "/Tables - Table 1"))
  write.csv(chem_subset, "chem_master_clean.csv", row.names = FALSE)
  
  #############################################################################################################
  #############################################################################################################
  #############################################################################################################
  
  # Reset the directory to the project folder
  setwd(current_directory)
}