#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
##########################################  CREATE THE LONG DATASET  ##########################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Purpose: This function creates and saves the long nhanes_subset_dataset that I'll be using for analysis
#          
# Inputs:   nhanes_subset              - dataframe containing all variables including cell types, chemicals,
#                                        and covariates (nhanes_subset_dataset)
#           subset_chemicals           - dataframe of chemicals to include for analysis
#
# Outputs:  long_nhanes_subset_dataset - dataframe of subsetted covariates, cell types, and chemicals

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

long_nhanes_subset_scale_function <- function(nhanes_subset,
                                              subset_chemicals)
{
  library(tidyverse)
  
  # setwd(current_directory)
  
  chems <- subset_chemicals$chemical_codename_use
  
  #############################################################################################################
  ################################ Create The Long Dataset For Chemicals And Cells ############################
  #############################################################################################################
 
  #make the long dataset - grouped by chemicals
  long_chemicals <- pivot_longer(data = nhanes_subset,
                                 cols = all_of(chems),
                                 names_to = "chemical_codename",
                                 values_to = "chem_measurement") %>%
    drop_na(chem_measurement)
  
  # print(str(long_chemicals))
  # print(unique(long_chemicals$chemical_codename))
  # print(length(unique(long_chemicals$SEQN)))

  #############################################################################################################

  immune <- c(
              "LBDLYMNO", #lymphocyte counts
              "LBDMONO",  #monocyte counts
              "LBDNENO",  #neutrophil counts
              "LBDEONO",  #eosinophil counts
              "LBDBANO",  #basophil counts
              "LBXWBCSI", #WBC count
              "LBXRBCSI", #RBC count
              "LBXMCVSI"  #MCV
  )
  
  #make the long dataset - grouped by chemicals and cell types
  long_chemicals_cells <- pivot_longer(data = long_chemicals,
                                       cols = all_of(immune),
                                       names_to = "celltype_codename",
                                       values_to = "cell_measurement")

  # print(str(long_chemicals_cells))

  rm(long_chemicals)

  #############################################################################################################
  ###################################### Create Log Transformed Variables #####################################
  #############################################################################################################

  #chem_log_measurement - also log transform creatinine and smoking variable
  long_nhanes_subset_dataset <- long_chemicals_cells %>%
    group_by(chemical_codename) %>%
    mutate(chem_log_measurement = scale(log2(chem_measurement))) %>%
    mutate(SMOKING = log2(SMOKING)) %>%
    mutate(URXUCR = log2(URXUCR)) %>%
    drop_na(chem_log_measurement) %>%
    dplyr::select(-chem_measurement) %>%
    mutate(RIAGENDR = relevel(factor(RIAGENDR), ref = "1"),
           RIDRETH1 = relevel(factor(RIDRETH1), ref = "3"))
  
  long_nhanes_subset_dataset$SDDSRVYR <- as.integer(long_nhanes_subset_dataset$SDDSRVYR)

  # print(str(long_nhanes_subset_dataset))

  #############################################################################################################
  return(long_nhanes_subset_dataset)
}