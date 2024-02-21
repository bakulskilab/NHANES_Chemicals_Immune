#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
########################################  CALCULATE SUMMARY STATISTICS  #######################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Purpose: This function calculates summary statistics for chemicals and cell-types
#          
# Inputs:   nhanes_subset      - dataframe of complete demographics, cells, and chemicals
# 
#           subset_chemicals   - dataframe of percent measurements above LOD per included chemical
#
# Outputs:  summary statistics - csv file containing all summary statistics for chemicals

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

table_chem_cell_stats <- function(nhanes_subset,
                                  subset_chemicals)
{
  library(tidyverse)
  
  #set the number of digits after the decimal point to be 2
  # options(digits=2)
  
  #TEMPORARY
  # nhanes_subset <- nhanes_subset_dataset
  # subset_chemicals <- use_these_chems

  setwd(current_directory)
  
  #############################################################################################################
  ############################################# Make Long Datasets ############################################
  #############################################################################################################
  
  chems <- subset_chemicals$chemical_codename_use
  
  #make the long dataset - grouped by chemicals, not log transformed
  long_chemicals <- gather(data = nhanes_subset, #this is the wide dataset subsetted to full demog and chems
                           key = chemical_codename_use, #this is the name of the new column to describe the chemicals
                           value = measurement, #these are the chemical measurements
                           all_of(chems), #these are the columns to adjust
                           factor_key=TRUE #keeps the columns in order
                          ) %>%
                    drop_na(measurement) #have to specify measurement because R can't look through everything
  
  print("long chemicals dataset")
  print(str(long_chemicals)) #full dataset: 1886172 obs. of  21 variables
  
  #############################################################################################################
  ######################################## Get Statistics on Chemicals ########################################
  #############################################################################################################
  
  #calculate stats on chemicals
  stats <- long_chemicals %>%
            group_by(chemical_codename_use) %>%
            summarise(
              n = n(), #total number of measurements
              min = round(min(measurement), digits = 2),
              quantile_25 = round(quantile(measurement, probs = 0.25), digits = 2),
              quantile_75 = round(quantile(measurement, probs = 0.75), digits = 2),
              max = round(max(measurement), digits = 2),
              mean = round(mean(measurement), digits = 2),
              stdev = round(sd(measurement), digits = 2),
              median = round(median(measurement), digits = 2)) %>%
            arrange(n) %>%
            ungroup()
  # print(str(stats))
  
  #add the chemical name to stats to make it understandable
  stats_chems <- left_join(subset_chemicals, stats, by = "chemical_codename_use")

  #reorder the table and drop some columns
  stats_chems <- stats_chems %>%
    dplyr::select(chemical_name,
                  chemical_codename_use,
                  n,
                  min,
                  quantile_25,
                  median,
                  quantile_75,
                  max,
                  mean,
                  stdev,
                  percent_above_LOD) %>%
    mutate(percent_above_LOD = round(percent_above_LOD, digits = 3))
  print("stats_chems shortened")
  print(dim(stats_chems))
  
  #save the table as a csv
  setwd(paste0(current_directory, "/Tables - Table 1"))
  write.csv(stats_chems, file = "subset_chemical_basic_statistics.csv", row.names = FALSE)
  setwd(current_directory)
  print("stats on chemicals saved as csv")

  #############################################################################################################
}