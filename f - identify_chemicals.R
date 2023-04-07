#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
##############################  Determine Which Chemicals To Include In Analysis  #############################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Purpose: This function selects which chemicals to use based on the limits of detection (LOD) per chemical
#          
# Inputs:   nhanes_full_dataset   - dataframe containing all variables including cell types, chemicals,
#                                   and covariates (nhanes_merged_dataset)
#           nhanes_comments       - dataframe of LOD comments for each chemical measurement (comments_clean)
#           chemical_dataset      - dataframe of chemical measurements
#           chem_master           - dataframe of information about the chemicals
#           weights_dataset       - dataframe of nhanes survey weights
#
# Outputs:  use_these_chems       - dataframe of chemicals with 50% of measurements above LOD

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

identify_chemicals <- function(nhanes_full_dataset,
                               nhanes_comments,
                               chemical_dataset,
                               chem_master,
                               weights_dataset)
{
  library(tidyverse)
  library(sjlabelled)
  library(magrittr) # %T>% package
  library(usefun)
  
  setwd(current_directory)
  
  #TEMPORARY
  # nhanes_full_dataset <- nhanes_merged_dataset
  # nhanes_comments <- comments_clean
  # chemical_dataset <- chemicals_clean
  # chem_master <- list_master_files$Chemicals
  # weights_dataset <- survey_weights

  #############################################################################################################
  
  #OUTLINE:
  # 1) Identify variables to start with (a. chemicals, b. immune measures, c. demographics)
  # 2) Subset participants dataset to only includes chemicals, immune measures, and demographics
  # 3) Drop urinary cadmium / molybdenum cross-reactive measurements, Drop NNAL measurement with creatinine = 0
  # 4) Drop participant measurements with weights of 0 - NOT NECESSARY
  # 5) Drop participants with no chemical measures
  # 6) Drop participants with urinary chemical measurements but no creatinine measurements
  # 7) Drop lipid unadjusted chemicals
  # 8) In comments dataset, keep only the participants with nonmissing data
  # 9) Drop chemicals with <50% of measurements above the LOD
  # 10) Clean up any issues

  #############################################################################################################
  ############################################# Identify Variables ############################################
  #############################################################################################################
  
  # Step 1: identify variables for a) chemicals, b) immune measures, and c) demographics
  
  # 1a) Chemicals
  # Make a vector of all chemical codenames, and drop SEQN and survey year
  total_chems <- chemical_dataset %>%
    dplyr::select(-SEQN,
                  -SDDSRVYR) %>%
    colnames(.)
  #537 chemicals (including creatinine, not including iron)
  
  # Identify and remove dietary components - not interested (34),
                      # phytoestrogens - also dietary (6),
                      # estradiol - not environmental (1),
                      # melamine - missing comments (1),
                      # cyanuric acid - missing comments (1),
                      # smoking chems - only measured in smokers (7),
                      # creatinine - I'll add it back with the demographics (1)
  chems <- chem_master %>%
    filter(!chem_family %in% c("Dietary Components",
                               "Phytoestrogens")) %>%
    filter(!chemical_codename_use %in% c("LBXEST",
                                          "SSCYA",
                                          "SSMEL",
                                          "URXUCR",
                                          "URXANBT",
                                          "URXANTT",
                                          "URXCOXT",
                                          "URXHPBT",
                                          "URXNICT",
                                          "URXNNCT",
                                          "URXNOXT")) %>%
    distinct(chemical_codename_use) %>%
    pull(chemical_codename_use)
  #486
  
  #############################################################################################################
  
  #1b/1c) Immune measures and Demographics
  
  #make a vector of the 8 immune measures
  immune <- c("LBXLYPCT", #lymphocytes
              "LBXMOPCT", #monocytes
              "LBXNEPCT", #neutrophils
              "LBXEOPCT", #eosinophils
              "LBXBAPCT", #basophils
              "LBXWBCSI", #WBC count
              "LBXRBCSI", #RBC count
              "LBXMCVSI")  #MCV
  
  #make a vector of the 5 demographics + SEQN + creatinine
  demog <- c("SEQN",     #ID
             "RIDAGEYR", #age
             "RIDRETH1", #race
             "RIAGENDR", #gender
             "INDFMPIR", #poverty-income ratio
             "BMXWAIST", #waist circumference
             "URXUCR")   #creatinine
  
  #############################################################################################################
  ########################################### Create Initial Dataset ##########################################
  #############################################################################################################
  
  # Step 2. Subset participants dataset to only includes chemicals, immune measures, and demographics
  
  #subset nhanes_full_dataset down to only these variables
  #exclude children (<18 years)
  #exclude participants with missing data
  #exclude participants with missing cotinine because it will become the smoking covariate
  
  nhanes_subset_unclean <- nhanes_full_dataset %>%
    dplyr::select(all_of(demog),
                  all_of(immune),
                  all_of(chems)) %T>%
    {print("total dataset")} %T>%
    {print(dim(.))} %>% #101316    501 (7 demog + 8 immune + 486 chems)
    
    #drop participants under 18 (children)
    filter(RIDAGEYR >= 18)  %T>%
    {print("drop <18 years")} %T>%
    {print(dim(.))} %>% #59204
    
    #dropping participants with missing demographics
    drop_na(RIDAGEYR,
            RIDRETH1,
            RIAGENDR) %T>%
    {print("drop age, race, gender")} %T>%
    {print(dim(.))} %>% #59204
    drop_na(INDFMPIR) %T>%
    {print("drop poverty-income ratio missings")} %T>%
    {print(dim(.))} %>% #53446
    drop_na(BMXWAIST) %T>%
    {print("drop waist circumference missings")} %T>%
    {print(dim(.))} %>% #48484
    drop_na(LBXCOT) %T>%
    {print("drop cotinine missings")} %T>%
    {print(dim(.))} %>% #46080
    
    #dropping participants with missing cell counts
    drop_na(LBXWBCSI) %T>%
    {print("drop WBC missings")} %T>%
    {print(dim(.))} %>% #46004
    drop_na(LBXRBCSI) %T>%
    {print("drop RBC missings")} %T>%
    {print(dim(.))} %>% #46004
    drop_na(LBXLYPCT,
            LBXMOPCT,
            LBXNEPCT,
            LBXEOPCT,
            LBXBAPCT,
            LBXMCVSI) %>%
    filter(!SEQN == "102389")  %T>% #WBC count of 400 which is > upper LOD
    {print("drop cell type missings")} %T>%
    {print(dim(.))} #45870
  
  print("dataset pre measurement cleaning")
  print(dim(nhanes_subset_unclean))
  #45870   501
  #this creates a dataset of all variables and all complete participants (good)
  
  #############################################################################################################
  ################################### Drop Urinary Cadmium Measurements of 0 ##################################
  #############################################################################################################
  
  # Step 3. Drop participants with urinary cadmium measurements of 0
  # Measurements are incorrect in NHANES due to urinary cadmium measurements being crossed with molybdenum
  
  #how many measurements
  sum(!is.na(nhanes_subset_unclean$URXUCD)) #13592 total measurements
  nhanes_subset_unclean %>%
    dplyr::select(URXUCD) %>%
    filter(URXUCD == 0) %>%
    nrow() #7
  
  #figure out who these participants are for dropping their measurements later in the comments dataset
  ucd_partic <- nhanes_subset_unclean %>%
    filter(URXUCD == 0) %>%
    pull(SEQN)
  
  # Convert measurements of 0 to NA
  nhanes_subset_cleaning <- nhanes_subset_unclean %>%
    mutate(URXUCD = ifelse(URXUCD == 0, NA, URXUCD))
  # 45870   501 (486 chems + 15 other variables)
  
  #check new measurement counts
  sum(!is.na(nhanes_subset_cleaning$URXUCD)) #13585 total measurements
  
  #############################################################################################################
  ################################# Drop Participants With Survey Weights of 0 ################################
  #############################################################################################################
  
  # Step 4. Drop participants with survey weight = 0
  #some participants have survey weights of 0 because the were "non-respondents" to a survey
  
  # #identify participants to initially keep
  # initial_seqn <- nhanes_subset_cleaning %>% pull(SEQN)
  # 
  # #keep only those participants in the weights dataset
  # weights_seqn_subset <- weights_dataset %>%
  #   filter(SEQN %in% initial_seqn)
  # identical(nhanes_subset_cleaning$SEQN, weights_seqn_subset$SEQN)
  # 
  # #how many 0s per chemical?
  # long_weights_temp <- weights_seqn_subset %>%
  #   pivot_longer(cols = WT_VNSUMDEHP:WT_URXIPM3,
  #                names_to = "chem_weights",
  #                values_to = "weights") %>%
  #   filter(weights == 0)
  # #no chemicals have weights of 0 so move on
  # 
  # rm(weights_seqn_subset, long_weights_temp)
  
  #Vy already converted the 0s to NAs when cleaning the weights dataset
  
  #############################################################################################################
  ############################## Drop Participants With No Chemical Measurements ##############################
  #############################################################################################################
  
  # Step 5. Drop participants with no chemical measurements
  
  print("______________")
  
  # Make a long dataset and drop participants
  long_cleaning <- nhanes_subset_cleaning %>%
    pivot_longer(cols = all_of(chems),
                 names_to = "chemical_codename",
                 values_to = "measurements") %>%
    drop_na(measurements) #%T>%
    # {print("drop measurements")} %T>%
    # {print(length(unique(.$SEQN)))} #no participants dropped
  
  #############################################################################################################
  ############################## Drop Participants With No Creatinine Measurements ############################
  #############################################################################################################
  
  # Step 6. Drop participants with urinary chemical measurements but no creatinine measurements
  
  long_creatinine_seqn <- long_cleaning %>%
    mutate(urinary = case_when(str_detect(chemical_codename, "^LB", negate = TRUE) ~ 99999999)) %>%
    mutate(creat_na = ifelse(urinary == 99999999 & is.na(URXUCR), NA, measurements)) %>%
    drop_na(creat_na) %>%
    dplyr::select(-urinary,
                  -creat_na) %>%
    distinct(SEQN) %>%
    pull(SEQN)
  
  rm(long_cleaning)
  length(long_creatinine_seqn) #45528 - dropped 342 participants
  
  #keep this subset of participants
  nhanes_creatinine <- nhanes_subset_cleaning %>%
    filter(SEQN %in% long_creatinine_seqn)
  
  
  
  
  #also drop URXNAL measurement from SEQN: 65564 because their creatinine measurement is log2(0) = -inf
  #Drop the creatinine=0 too because this participant has no urinary measures in this subset
  #Only 1 person has creatinine=0
  nhanes_creatinine <- nhanes_creatinine %>%
    mutate(URXNAL = case_when(SEQN == 65564 & URXNAL == 0.0553 ~ 5000000000,
                              TRUE ~ as.numeric(URXNAL))) %>%
    mutate(URXNAL = na_if(URXNAL, 5000000000)) %>%
    mutate(URXUCR = ifelse(URXUCR == 0, NA, URXUCR))
    
    
  #############################################################################################################
  ###################################### Drop Lipid Unadjusted Chemicals ######################################
  #############################################################################################################
  
  # Step 7. Drop lipid unadjusted chemicals
  
  #select the lipid adjusted chemicals (n = 79)
  chem_lipids <- chem_master %>%
    filter(grepl("adj", chemical_name) | 
             grepl("Adj", chemical_name) |
             grepl("lipid ad", chemical_name)) %>%
    distinct(chemical_codename_use, .keep_all = TRUE)
  
  #remove the LA ending from the codenames to identify the unadjusted codenames
  unadjusted_lipid_codename <- gsub("LA$", "", chem_lipids$chemical_codename_use)
  #rename the column to chemical_codename for merging purposes
  unique_unadjusted_lipid <- unique(unadjusted_lipid_codename)
  #79 lipid-unadjusted chemicals
  
  #PCB 199 codename is just different than the rest of the PCBs
  unique_unadjusted_lipid <- gsub("LBX199", "LBD199", unique_unadjusted_lipid)
  
  #drop the unadjusted duplicates of chemicals
  nhanes_subset_lipid <- nhanes_creatinine %>%
    dplyr::select(-all_of(unique_unadjusted_lipid))
  print("dim nhanes subset with only lipids")
  print(dim(nhanes_subset_lipid))
  #45528   422 (407 chems + 15 other variables)
  
  #############################################################################################################
  ########################################## Subset Comments Dataset ##########################################
  #############################################################################################################
  
  # Step 8. In comments dataset, keep only the participants with nonmissing data
  
  #in the comments dataset, keep only the participants I'm going to include
  comments_participants <- nhanes_comments %>%
    filter(SEQN %in% long_creatinine_seqn) %>%
    dplyr::select(-SDDSRVYR)
  # print(dim(comments_participants))
  #45528   427
  
  
  #make vector of comment_codenames - excluding SEQN
  num_columns <- ncol(comments_participants)
  comment_codenames <- colnames(comments_participants)[2:num_columns]
  
  
  #make the comments dataset into long form by chemical
  #also drop those urinary cadmium-molybdenum contaminated measurements
  #also drop the nnal participant who has a creatinine measurement of 0
  long_comments_sub <- pivot_longer(data = comments_participants,
                                    cols = all_of(comment_codenames),
                                    names_to = "comment_codename",
                                    values_to = "comment") %>%
    mutate(comment = ifelse(SEQN %in% ucd_partic & comment_codename == "URDUCDLC", NA, comment)) %>%
    mutate(comment = ifelse(SEQN == 65564 & comment_codename == "URDNALLC", NA, comment)) %>%
    dplyr::select(-SEQN)
  
  # print(table(long_comments_sub$comment))
  
  #############################################################################################################
  ######################################## Exclude Chemicals Under LOD ########################################
  #############################################################################################################
  
  # Step 9. Drop chemicals with <50% of measurements above the LOD
  
  counts_stats <- long_comments_sub %>%
    group_by(comment_codename) %>%
    dplyr::summarise(above = sum(comment == 0 | comment == 2, na.rm = TRUE),
                     below = sum(comment == 1, na.rm = TRUE),
                     total = above + below,
                     percent_above_LOD = (above / total)*100) %>%
    mutate(included =
             ifelse(percent_above_LOD >= 50, "yes", "no")) %>%
    ungroup()
  
  
  # print(dim(counts_stats)) #425   6
  
  sum(counts_stats$percent_above_LOD >= 50)
  #215 (7 smk + 6 phyto + 1 estradiol + 5 dietary + 196 chemicals)
  
  
  #clean up the chemical info dataset
  weird_smk_chems <- c("URXANBT",
                       "URXANTT",
                       "URXCOXT",
                       "URXHPBT",
                       "URXNICT",
                       "URXNNCT",
                       "URXNOXT")
  
  #select only the useful columns
  chem_conv <- chem_master %>%
    dplyr::select(chemical_codename_use,
                  chemical_name,
                  comments_codename_use,
                  chem_family,
                  chem_family_shortened) %>%
    distinct(chemical_codename_use, .keep_all = TRUE) %>%
    filter(!chem_family == "Dietary Components") %>%
    #creatinine doesn't have comments
    filter(!chemical_codename_use == "URXUCR") %>%
    #estradiol and phytoestrogens are not environmental exposures
    filter(!chemical_codename_use == "LBXEST",
           !chem_family == "Phytoestrogens") %>%
    #melamine and cyanuric acid don't have comments available
    filter(!chemical_codename_use == "SSMEL",
           !chemical_codename_use == "SSCYA") %>%
    #drop the weird smoking chemicals
    filter(!chemical_codename_use %in% weird_smk_chems) %>%
    #drop unadjusted lipid chemicals
    filter(!chemical_codename_use %in% unique_unadjusted_lipid) %>%
    dplyr::rename(comment_codename = comments_codename_use)
  #407 chemicals (no lipids, phytoestrogens, estradiol, melamine, cyanuric acid, those 7 smoking problem chems, or creatinine)
  
  #merge counts_stats with chemicals
  chem_comms <- left_join(chem_conv, counts_stats, by = "comment_codename")
  print(dim(chem_comms))
  #407 10
  
  #drop chemicals with 50% measurements below LOD
  chem_comms_above <- chem_comms %>%
    filter(percent_above_LOD >= 50)
  #196 - missing copper and zinc
  
  #############################################################################################################
  ############################################ Clean Up The Dataset ###########################################
  #############################################################################################################
  
  #get measurement counts for copper
  total_counts <- length(nhanes_subset_lipid$LBXSCU)
  miss_cu <- sum(is.na(nhanes_subset_lipid$LBXSCU))
  total_cu <- total_counts - miss_cu
  
  #get measurement counts for zinc
  total_counts <- length(nhanes_subset_lipid$LBXSZN)
  miss_zn <- sum(is.na(nhanes_subset_lipid$LBXSZN))
  total_zn <- total_counts - miss_zn
  
  #manually add copper and zinc back into the dataset
  copper <- c("LBXSCU",
              "Serum Copper (ug/dL)",
              NA,
              "Metals",
              "Metals",
              total_cu, #above
              "0",      #below
              total_cu, #total
              "100",    #percent above
              NA)
  zinc <- c("LBXSZN",
            "Serum Zinc (ug/dL)",
            NA,
            "Metals",
            "Metals",
            total_zn, #above
            "0",      #below
            total_zn, #total
            "100",    #percent above
            NA)
  
  #turn the vectors into a dataframe
  copper_zinc <- as.data.frame(t(data.frame(copper, zinc)))
  #add column names
  chem_colnames <- colnames(chem_comms_above)
  colnames(copper_zinc) <- chem_colnames
  #remove rownames
  copper_zinc <- remove_rownames(copper_zinc)
  
  #rbind the dataframe onto the bottom of chem_comms_above
  chem_bind <- rbind(chem_comms_above, copper_zinc)
  
  #use these chemicals
  use_these_chems <- chem_bind %>%
    dplyr::select(chemical_codename_use,
                  chemical_name,
                  chem_family,
                  chem_family_shortened,
                  comment_codename,
                  above,
                  total,
                  percent_above_LOD)
  #198 chemicals for analysis
  
  #reset the data types
  use_these_chems$above <- as.integer(use_these_chems$above)
  use_these_chems$total <- as.integer(use_these_chems$total)
  use_these_chems$percent_above_LOD <- as.numeric(use_these_chems$percent_above_LOD)
  
  print(paste("number of chemicals with LOD > 50%:", length(use_these_chems$chemical_codename_use)))
  
  #############################################################################################################
  #############################################################################################################
  #############################################################################################################
  
  # Reset the directory to the project folder
  setwd(current_directory)
  
  return(use_these_chems)
}