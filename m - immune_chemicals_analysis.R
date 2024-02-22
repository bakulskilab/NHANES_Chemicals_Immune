#############################################################################################################
#######################   MAIN SCRIPT - ENVIRONMENTAL CHEMICALS AND IMMUNE MEASURES   #######################
#############################################################################################################

#set working directory
current_directory <- "C:/Users/lamid/OneDrive/Documents/Research/Chemicals-and-Cell-Proportions-NHANES"
# current_directory <- "C:/Users/T7920/Desktop/Lauren/NHANES Project"
setwd(current_directory)

#TEMPORARY
# comments_clean <- readRDS("comments_clean.rds")
# chemicals_clean <- readRDS("chemicals_clean.rds")
# demographics_clean <- readRDS("demographics_clean.rds")
# list_master_files <- readRDS("list_master_files.rds")
# response_clean <- readRDS("response_clean.rds")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  Download and Load Packages  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

source("f - download_packages.R", local=T)
download_packages()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  Set Up Folders  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

source("f - create_folders.R", local=T)
create_folders()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  Merge the Clean Datasets  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

chemicals_clean <- chemicals_clean %>%
  dplyr::select(-study_year)

source("f - merge_datasets_together.R", local=T)
nhanes_merged_dataset <- merge_datasets_together(demographics_dataset = demographics_clean
                                                 , response_dataset = response_clean
                                                 , chemicals_dataset = chemicals_clean)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Determine Which Chemicals To Include In Analysis ~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# source("f - limit_of_detection.R", local=T)
source("f - identify_chemicals.R", local=T)
use_these_chems <- identify_chemicals(nhanes_full_dataset = nhanes_merged_dataset,
                                      nhanes_comments = comments_clean,
                                      chemical_dataset = chemicals_clean,
                                      chem_master = list_master_files$Chemicals,
                                      weights_dataset = survey_weights)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Create The Working Dataset ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

source("f - nhanes_subset_function.R", local=T)
nhanes_subset_dataset <- nhanes_subset_function(nhanes_full_dataset = nhanes_merged_dataset,
                                                subset_chemicals = use_these_chems)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Create The Long Working Dataset ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#log2 transformed
source("f - long_nhanes_subset_function.R", local=T)
long_nhanes_subset_dataset <- long_nhanes_subset_function(nhanes_subset = nhanes_subset_dataset,
                                                          subset_chemicals = use_these_chems)
#log2 and standardized
source("f - long_nhanes_subset_scale_function.R", local=T)
long_nhanes_subset_scale_dataset <- long_nhanes_subset_scale_function(nhanes_subset = nhanes_subset_dataset,
                                                                      subset_chemicals = use_these_chems)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Create Barplot and Scatterplot of LOD data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#supplemental figures 1 and 2 - rerun 8/8/23, 2/21/24
source("f - limit_of_detection_bar_plots.R", local=T)
limit_of_detection_bar_plots(subset_chemicals = use_these_chems)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  Make Table 1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#tables 1 and 2 - unweighted
source("f - table_demographics.R", local=T)
table_demographics(nhanes_subset = nhanes_subset_dataset)

#table 1 - weighted - manually fix IQRs
source("f - table_1_weighted.R", local=T)
table_1_weighted(nhanes_subset = nhanes_subset_dataset,
                 demog_dataset = demographics_clean)

#table 2 - weighted - manually fix IQRs and add adult normal range
source("f - table_2_weighted.R", local=T)
table_2_weighted(nhanes_subset = nhanes_subset_dataset,
                 demog_dataset = demographics_clean)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~  Put Together Chemicals File For Supplemental Table ~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#supplemental table 1 - reran 8/8/23
source("f - supplemental_table_chemicals.R", local = TRUE)
supplemental_table_chemicals(chem_master = list_master_files$Chemicals,
                             conversion = use_these_chems)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Make Inclusion/Exclusion Table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#supplemental table 2
source("f - table_inclusion_exclusion.R", local=T)
table_inclusion_exclusion(nhanes_subset = nhanes_subset_dataset,
                          nhanes_full_dataset = nhanes_merged_dataset,
                          demog_dataset = demographics_clean)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  Make Supplemental Table Chemicals ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#supplemental table 3 - updated 8/8/23
source("f - table_chem_cell_stats.R", local=T)
table_chem_cell_stats(nhanes_subset = nhanes_subset_dataset,
                      subset_chemicals = use_these_chems)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  Make Supplemental Table Immune ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# table 2
source("f - table_immune_category_stats.R", local=T)
table_immune_category_stats(nhanes_subset = nhanes_subset_dataset)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Make Correlation Plots ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#figure 2 and supplemental table 4 - updated 8/8/23
source("f - correlation_plot_chemicals.R", local = TRUE)
correlation_plot_chemicals(subset_chemicals = use_these_chems,
                           nhanes_subset = nhanes_subset_dataset,
                           conversion = use_these_chems)

#supplemental figure 3 and supplemental table 5
source("f - correlation_plot_immune.R", local = TRUE)
correlation_plot_immune(nhanes_subset = nhanes_subset_dataset)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Calculate Correlation Stats ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#correlations of chemicals within chemical classes
source("f - correlation_stats.R", local = TRUE)
correlation_stats(subset_chemicals = use_these_chems,
                  nhanes_subset = nhanes_subset_dataset)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~  Weighted Linear Regression, Adjusted for Cotinine ~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#unscaled for written results
source("f - run_if_else_glm_weighted.R", local=T)
source("f - run_linear_regression.R", local=T)
model_stats_smk <- run_linear_regression(long_nhanes_subset = long_nhanes_subset_dataset,
                                         conversion = use_these_chems,
                                         weights_dataset = survey_weights,
                                         nhanes_subset = nhanes_subset_dataset)

setwd(paste0(current_directory, "/Regression Results"))
write.csv(model_stats_smk_clean, "model_stats_smk.csv", row.names = FALSE)
setwd(current_directory)

#scaled for figures
source("f - run_if_else_glm_weighted.R", local=T)
source("f - run_linear_regression.R", local=T)
model_stats_smk_scaled <- run_linear_regression(long_nhanes_subset = long_nhanes_subset_scale_dataset,
                                                conversion = use_these_chems,
                                                weights_dataset = survey_weights,
                                                nhanes_subset = nhanes_subset_dataset)

setwd(paste0(current_directory, "/Regression Results"))
write.csv(model_stats_smk_scaled, "model_stats_smk_scaled.csv", row.names = FALSE)
setwd(current_directory)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Interpret Beta Coefficients ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

source("f - interpret_beta.R", local=T)
interpret_beta(model_stats = model_stats_smk)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~  Weighted Linear Regression, Unadjusted for Cotinine ~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

source("f - run_if_else_glm_wt_nosmk.R", local=T)
source("f - run_linear_regression_no_smk.R", local=T)
model_stats_scale_no_smk <- run_linear_regression_no_smk(long_nhanes_subset = long_nhanes_subset_scale_dataset,
                                                         conversion = use_these_chems,
                                                         weights_dataset = survey_weights,
                                                         nhanes_subset = nhanes_subset_dataset)

setwd(paste0(current_directory, "/Regression Results"))
write.csv(model_stats_scale_no_smk, "model_stats_scale_no_smk.csv", row.names = FALSE)
setwd(current_directory)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Volcano Plot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#supp figure 4
source("f - volcano_plot_lin_reg.R", local = TRUE)
volcano_plot_lin_reg(model_stats = model_stats_smk_scaled,
                     conversion = use_these_chems)

#figure 3
source("f - significant_linear_chemical_plot.R", local = TRUE)
significant_linear_chemical_plot(model_stats = model_stats_smk_scaled,
                                 conversion = use_these_chems)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Forest Plots ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#figure 4 - reran 8/8/23
source("f - forest_plot_lin_reg.R", local = TRUE)
forest_plot_lin_reg(model_stats = model_stats_smk_scaled,
                    conversion = use_these_chems)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Make Correlation Plots ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#supplemental figure 5 - rerun 8/8/23
source("f - correlation_plot_unweight_weight.R", local = TRUE)
correlation_plot_unweight_weight(model_stats_wt_adj = model_stats_smk_scaled,
                                 model_stats_wt_unadj = model_stats_scale_no_smk)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Run Demographics Regressions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#linear regressions unscaled because no chemical covariates
source("f - run_demog_regression.R", local=T)
run_demog_regression(nhanes_subset = nhanes_subset_dataset,
                     conversion = use_these_chems,
                     demog_dataset = demographics_clean)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Run Logistic Regressions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#unscaled for written results
source("f - run_logit_regression.R", local=T)
run_logit_regression(long_nhanes_subset = long_nhanes_subset_dataset,
                     conversion = use_these_chems,
                     weights_dataset = survey_weights)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#