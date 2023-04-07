#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
##############################################  Set Up Folders  ###############################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Purpose: This function creates the folders for the whole project
#          
# Inputs:   NA
#
# Outputs:  Folders created in computer directory

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

create_folders <- function()
{
  # NOTE: Current directory should be updated and set in main script first
  current_directory <- getwd()
  
  # Determine all file names in the current working directory
  all_files_in_current_directory <- list.files()
  
  #Set a vector of folder names
  name_of_folder <- c("Bar Plots - Sample Size by Toxicant",
                      "Correlation Plots - Demog, Cells, Chemicals",
                      "Forest Plot - Chemicals by Cell Types",
                      "Regression Results",
                      "Tables - Table 1",
                      "Volcano Plots")
  
  
  # Make a new folder if the folder doesn't exist and print location
  for(x_folder in name_of_folder)
  {
    if(file.exists(x_folder))
    {
      print("Folder exists - do nothing")
      new_working_directory <- paste0(current_directory, "/", x_folder)
      print(new_working_directory)
    } else {
      print("Folder NOT exist - create folder")
      dir.create(x_folder)
      new_working_directory <- paste0(current_directory, "/", x_folder)
      print(new_working_directory)
    }
  }
  
  
}