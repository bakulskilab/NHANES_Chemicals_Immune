#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#########################################  DOWNLOAD AND LOAD PACKAGES  ########################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Purpose: This function installs and loads the packages
#          
# Inputs: None
#
# Outputs: loaded packages for this project

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

download_packages <- function()
{
  # Create vector of packages to install for this project
  packages_to_install <- c("tidyverse",
                           "usefun",
                           "sjlabelled",
                           "magrittr",
                           "ggplot2",
                           "gtsummary",
                           "compareGroups",
                           "flextable",
                           "dichromat",
                           "pheatmap",
                           "corrplot",
                           "broom",
                           "survey",
                           "ggrepel",
                           "rstatix",
                           "gt",
                           "colorspace",
                           "ggpubr")
  
  # Install the packages
  install.packages(packages_to_install)
  
  # Load the packages
  for (pkg in packages_to_install) {
    library(pkg, character.only = TRUE)
  }
  
  rm(pkg)
  
  library(devtools)
  install_github("jokergoo/ComplexHeatmap")
  
  print("packages loaded")
}