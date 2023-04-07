#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
################################  FUNCTION TO MERGE NHANES DATASET TOGETHER  ##################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Purpose: This function merge all NHANES dataset together and include only variables that are need
#          
# Inputs: demographics_dataset - data frame containing the demographic features for each participant 
#        # mortality_dataset <- data frame containing information on deceased status and time until death
#         response_dataset <- data frame containing measurements for 10 allostatic load components
#         chemicals_dataset <- data frame containing measurements for >400 chemical biomarkers
#         weights_dataset <- data frame containing the values on survey weights
#
# Outputs: dataset_merged_updated - clean merged dataset of demographics, mortality, response, and weights datasets

merge_datasets_together <- function(demographics_dataset
                                    # , weights_dataset
                                    # , mortality_dataset
                                    , chemicals_dataset
                                    , response_dataset)
{
  #load libraries
  library(tidyverse)
  library(sjlabelled)
  
  # Determine the number of arguments
  num_of_arguments <- nargs()
  
  # Drop SDDSRVYR from two datasets
  demographics_dataset <- demographics_dataset #%>%
    # dplyr::select(-SDDSRVYR)
  response_dataset <- response_dataset %>%
    dplyr::select(-SDDSRVYR)
  chemicals_dataset <- chemicals_dataset %>%
    dplyr::select(-SDDSRVYR)
  
  # Assemble a list containing the datasets
  list_of_datasets <- list(demographics_dataset = demographics_dataset
                           # , weights_dataset = weights_dataset
                           , response_dataset = response_dataset  
                           , chemicals_dataset = chemicals_dataset)
  
  # Define a function to join the dataset by the participant identifiers
  joining_by_seqn <- function(x, y) full_join(x, y, by = "SEQN")

  # Merge all the datasets together
  dataset_merged <- list_of_datasets %>%
    reduce(joining_by_seqn)


  # Copy the attributes and labels into the merged dataset
  for(i in seq(num_of_arguments))
  {
    dataset_merged <- copy_labels(df_new = dataset_merged, df_origin = list_of_datasets[[i]])
  }

  # Define a vector of the variable codename to include
  colnames_include <- colnames(dataset_merged)

  # Define a vector of codenames that are not duplicated
  colnames_include_updated <- colnames_include[grepl("\\.1$", colnames_include) == FALSE]

  # Define a dataframe of pertinent variables for further analyses
  dataset_merged_updated <- dataset_merged[,colnames_include_updated]

  return(dataset_merged_updated)
}