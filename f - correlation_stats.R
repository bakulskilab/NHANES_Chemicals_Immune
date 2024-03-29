#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
##############################  MAKE PLOT OF CORRELATIONS BETWEEN THE CHEMICALS  ##############################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Purpose: This function creates correlation plots of demographics, chemicals, and immune measures
#          
# Inputs: nhanes_subset     - dataframe containing complete demographic and immune measures data for each
#                             participant
#         conversion        - dataframe to convert between chemical_codename and chemical family
#         
#         subset_chemicals  - dataframe that has a vector of the chemicals to use
#
# Outputs: Plot showing correlations between chemicals

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

correlation_stats <- function(subset_chemicals,
                              nhanes_subset)
{
  library(tidyverse)
  # library(gplots) #heatmap.2 package
  # library(Hmisc)
  # library(data.table)
  
  # subset_chemicals <- use_these_chems
  # nhanes_subset <- nhanes_subset_dataset
  
  #############################################################################################################
  ######################### Clean The Chemical Data To Make It Usable For Correlations ########################
  #############################################################################################################
  
  #select the usable chemicals
  chems <- subset_chemicals$chemical_codename_use
  
  #select those chemicals from the wide nhanes dataset,
  #log2 transform values
  nhanes_subset_chems <- nhanes_subset %>%
    dplyr::select(all_of(chems)) %>%
    log2(.)
  
  #############################################################################################################
  ####################################### Calculate Chemical Correlations #####################################
  #############################################################################################################
  
  #calculate the correlations
  chem_correlations <- cor(nhanes_subset_chems,
                           use = "pairwise.complete.obs",
                           method = "spearman"
                           # ,
                           # diag = FALSE
                           )
  
  # identical(colnames(chem_correlations), colnames(nhanes_subset_chems))
  #TRUE - the correlations are in the same order as the chem families dataset


  
  # info from: https://sites.ualberta.ca/~ahamann/teaching/graphics/LabHM.pdf
  
  #since the correlations are outputted as a duplicated square matrix rather than triangle, 
  #need to select only one triangle of data
  lower_tri <- function(chem_correlations){chem_correlations[upper.tri(chem_correlations)] <- NA
  return(chem_correlations)}
  bottom_tri_correlations <- as.data.frame(lower_tri(chem_correlations))
  
  #remove diagonal perfect correlations
  diag(bottom_tri_correlations) <- NA
  #############################################################################################################
  ##################################### Calculate Chem Family Correlations ####################################
  #############################################################################################################
  
  cor_family <- data.frame(matrix(ncol = 17, nrow = 1,
                                  dimnames = list(NULL)))
  
  chem_family_levels <- c("Acrylamide"
                          # , "Melamine"
                          , "Brominated Flame Retardants (BFR)"
                          , "Phosphate Flame Retardants (PFR)"
                          , "Polychlorinated Biphenyls (PCB)"
                          , "Dioxins"
                          , "Furans"
                          , "Metals"
                          , "Phthalates & Plasticizers"
                          , "Personal Care & Consumer Product Compounds"
                          , "Pesticides"
                          , "Aromatic Amines"
                          # , "Phytoestrogens"
                          , "Polyaromatic Hydrocarbons (PAH)"
                          , "Volatile Organic Compounds (VOC)"
                          , "Smoking Related Compounds"
                          , "Per- and Polyfluoroalkyl Substances (PFAS)"
                          , "Aldehydes"
                          # , "Dietary Components"
                          , "Other")
  
  colnames(cor_family) <- chem_family_levels

  
  #acrylamide
  acr_codenames <- subset_chemicals %>%
    filter(chem_family == "Acrylamide") %>%
    pull(chemical_codename_use)
  acr_cor <- bottom_tri_correlations %>%
    dplyr::select(all_of(acr_codenames)) %>%
    rownames_to_column(var = "chem_rownames") %>%
    subset(., chem_rownames %in% acr_codenames) %>%
    dplyr::select(-chem_rownames)
  print("Acrylamide correlation")
  cor_family$Acrylamide <- mean(as.matrix(acr_cor), na.rm = TRUE)
  
  #BFRs - yes, using pull would have been easier, but I'm not fixing the whole thing
  bfr <- subset_chemicals %>%
    filter(chem_family == "Brominated Flame Retardants (BFR)")
  bfr_codenames <- bfr$chemical_codename_use
  bfr_cor <- bottom_tri_correlations %>%
    dplyr::select(all_of(bfr_codenames)) %>%
    rownames_to_column(var = "chem_rownames") %>%
    subset(., chem_rownames %in% bfr_codenames) %>%
    dplyr::select(-chem_rownames)
  print("BFR correlation")
  cor_family$`Brominated Flame Retardants (BFR)` <- mean(as.matrix(bfr_cor), na.rm = TRUE)
  
  #PFRs
  pfr <- subset_chemicals %>%
    filter(chem_family == "Phosphate Flame Retardants (PFR)")
  pfr_codenames <- pfr$chemical_codename_use
  pfr_cor <- bottom_tri_correlations %>%
    dplyr::select(all_of(pfr_codenames)) %>%
    rownames_to_column(var = "chem_rownames") %>%
    subset(., chem_rownames %in% pfr_codenames) %>%
    dplyr::select(-chem_rownames)
  print("PFR correlation")
  cor_family$`Phosphate Flame Retardants (PFR)` <- mean(as.matrix(pfr_cor), na.rm = TRUE)
  
  #PCBs
  pcb <- subset_chemicals %>%
    filter(chem_family_shortened == "PCBs")
  pcb_codenames <- pcb$chemical_codename_use
  pcb_cor <- bottom_tri_correlations %>%
    dplyr::select(all_of(pcb_codenames)) %>%
    rownames_to_column(var = "chem_rownames") %>%
    subset(., chem_rownames %in% pcb_codenames) %>%
    dplyr::select(-chem_rownames)
  print("PCB correlation")
  cor_family$`Polychlorinated Biphenyls (PCB)` <- mean(as.matrix(pcb_cor), na.rm = TRUE)
  
  #Dioxins
  dio <- subset_chemicals %>%
    filter(chem_family == "Dioxins")
  dio_codenames <- pfr$chemical_codename_use
  dio_cor <- bottom_tri_correlations %>%
    dplyr::select(all_of(dio_codenames)) %>%
    rownames_to_column(var = "chem_rownames") %>%
    subset(., chem_rownames %in% dio_codenames) %>%
    dplyr::select(-chem_rownames)
  print("Dioxin correlation")
  cor_family$`Dioxins` <- mean(as.matrix(dio_cor), na.rm = TRUE)
  
  #Furans
  furan <- subset_chemicals %>%
    filter(chem_family == "Dioxins")
  fur_codenames <- furan$chemical_codename_use
  fur_cor <- bottom_tri_correlations %>%
    dplyr::select(all_of(fur_codenames)) %>%
    rownames_to_column(var = "chem_rownames") %>%
    subset(., chem_rownames %in% fur_codenames) %>%
    dplyr::select(-chem_rownames)
  print("Furan correlation")
  cor_family$`Furans` <- mean(as.matrix(fur_cor), na.rm = TRUE)
  
  #Metals
  metal <- subset_chemicals %>%
    filter(chem_family == "Metals")
  metal_codenames <- metal$chemical_codename_use
  metal_cor <- bottom_tri_correlations %>%
    dplyr::select(all_of(metal_codenames)) %>%
    rownames_to_column(var = "chem_rownames") %>%
    subset(., chem_rownames %in% metal_codenames) %>%
    dplyr::select(-chem_rownames)
  print("Metals correlation")
  cor_family$`Metals` <- mean(as.matrix(metal_cor), na.rm = TRUE)
  print(cor_family$`Metals`)
  
  #P&P
  pp <- subset_chemicals %>%
    filter(chem_family == "Phthalates & Plasticizers")
  pp_codenames <- pp$chemical_codename_use
  pp_cor <- bottom_tri_correlations %>%
    dplyr::select(all_of(pp_codenames)) %>%
    rownames_to_column(var = "chem_rownames") %>%
    subset(., chem_rownames %in% pp_codenames) %>%
    dplyr::select(-chem_rownames)
  print("P&P correlation")
  cor_family$`Phthalates & Plasticizers` <- mean(as.matrix(pp_cor), na.rm = TRUE)
  
  #products
  prod <- subset_chemicals %>%
    filter(chem_family == "Personal Care & Consumer Product Compounds")
  prod_codenames <- prod$chemical_codename_use
  prod_cor <- bottom_tri_correlations %>%
    dplyr::select(all_of(prod_codenames)) %>%
    rownames_to_column(var = "chem_rownames") %>%
    subset(., chem_rownames %in% prod_codenames) %>%
    dplyr::select(-chem_rownames)
  print("prod correlation")
  cor_family$`Personal Care & Consumer Product Compounds` <- mean(as.matrix(prod_cor), na.rm = TRUE)
  
  #Pesticides
  pest <- subset_chemicals %>%
    filter(chem_family == "Pesticides")
  pest_codenames <- pest$chemical_codename_use
  pest_cor <- bottom_tri_correlations %>%
    dplyr::select(all_of(pest_codenames)) %>%
    rownames_to_column(var = "chem_rownames") %>%
    subset(., chem_rownames %in% pest_codenames) %>%
    dplyr::select(-chem_rownames)
  print("Pesticides correlation")
  cor_family$`Pesticides` <- mean(as.matrix(pest_cor), na.rm = TRUE)
  
  #Aromatic amines
  aro <- subset_chemicals %>%
    filter(chem_family == "Aromatic Amines")
  aro_codenames <- aro$chemical_codename_use
  aro_cor <- bottom_tri_correlations %>%
    dplyr::select(all_of(aro_codenames)) %>%
    rownames_to_column(var = "chem_rownames") %>%
    subset(., chem_rownames %in% aro_codenames) %>%
    dplyr::select(-chem_rownames)
  print("Aro correlation")
  cor_family$`Aromatic Amines` <- mean(as.matrix(aro_cor), na.rm = TRUE)
  
  #PAHs
  pah <- subset_chemicals %>%
    filter(chem_family == "Polyaromatic Hydrocarbons (PAH)")
  pah_codenames <- pah$chemical_codename_use
  pah_cor <- bottom_tri_correlations %>%
    dplyr::select(all_of(pah_codenames)) %>%
    rownames_to_column(var = "chem_rownames") %>%
    subset(., chem_rownames %in% pah_codenames) %>%
    dplyr::select(-chem_rownames)
  print("PAH correlation")
  print(mean(as.matrix(pah_cor), na.rm = TRUE))
  cor_family$`Polyaromatic Hydrocarbons (PAH)` <- mean(as.matrix(pah_cor), na.rm = TRUE)
  
  #VOCs
  voc <- subset_chemicals %>%
    filter(chem_family == "Volatile Organic Compounds (VOC)")
  voc_codenames <- voc$chemical_codename_use
  voc_cor <- bottom_tri_correlations %>%
    dplyr::select(all_of(voc_codenames)) %>%
    rownames_to_column(var = "chem_rownames") %>%
    subset(., chem_rownames %in% voc_codenames) %>%
    dplyr::select(-chem_rownames)
  print("VOC correlation")
  print(mean(as.matrix(voc_cor), na.rm = TRUE))
  cor_family$`Volatile Organic Compounds (VOC)` <- mean(as.matrix(voc_cor), na.rm = TRUE)
  
  #Smoking related compounds
  smoke <- subset_chemicals %>%
    filter(chem_family == "Smoking Related Compounds")
  smoke_codenames <- smoke$chemical_codename_use
  smoke_cor <- bottom_tri_correlations %>%
    dplyr::select(all_of(smoke_codenames)) %>%
    rownames_to_column(var = "chem_rownames") %>%
    subset(., chem_rownames %in% smoke_codenames) %>%
    dplyr::select(-chem_rownames)
  print("Smoking correlation")
  print(mean(as.matrix(smoke_cor), na.rm = TRUE))
  cor_family$`Smoking Related Compounds` <- mean(as.matrix(smoke_cor), na.rm = TRUE)

  #pfas
  pfas <- subset_chemicals %>%
    filter(chem_family == "Per- and Polyfluoroalkyl Substances (PFAS)")
  pfas_codenames <- pfas$chemical_codename_use
  pfas_cor <- bottom_tri_correlations %>%
    dplyr::select(all_of(pfas_codenames)) %>%
    rownames_to_column(var = "chem_rownames") %>%
    subset(., chem_rownames %in% pfas_codenames) %>%
    dplyr::select(-chem_rownames)
  print("pfas correlation")
  cor_family$`Per- and Polyfluoroalkyl Substances (PFAS)` <- mean(as.matrix(pfas_cor), na.rm = TRUE)
 
  #aldehydes
  ald <- subset_chemicals %>%
    filter(chem_family == "Aldehydes")
  ald_codenames <- ald$chemical_codename_use
  ald_cor <- bottom_tri_correlations %>%
    dplyr::select(all_of(ald_codenames)) %>%
    rownames_to_column(var = "chem_rownames") %>%
    subset(., chem_rownames %in% ald_codenames) %>%
    dplyr::select(-chem_rownames)
  print("aldehydes correlation")
  cor_family$`Aldehydes` <- mean(as.matrix(ald_cor), na.rm = TRUE)
  
  #other
  other <- subset_chemicals %>%
    filter(chem_family == "Other")
  other_codenames <- other$chemical_codename_use
  other_cor <- bottom_tri_correlations %>%
    dplyr::select(all_of(other_codenames)) %>%
    rownames_to_column(var = "chem_rownames") %>%
    subset(., chem_rownames %in% other_codenames) %>%
    dplyr::select(-chem_rownames)
  print("other correlation")
  cor_family$`Other` <- mean(as.matrix(other_cor), na.rm = TRUE)
  
  View(cor_family)
  
  #save the table as a csv
  setwd(paste0(current_directory, "/Correlation Plots - Demog, Cells, Chemicals"))
  write.csv(cor_family, file = "chemical_family_correlations.csv", row.names = FALSE)
  setwd(current_directory)
  print("correlation statistics saved as csv in Correlation Plots folder")
  
  #mean of family correlations
  print("mean of chemical family mean correlations")
  print(rowMeans(cor_family))
  #0.450584
  
  cor_family_t <- as.data.frame(t(cor_family)) %>% rename(correlations = V1)
  #11/17 families have correlation <0.5
  
  print("mean overall (without perfect correlations)")
  diag(bottom_tri_correlations) <- NA
  mean(unlist(bottom_tri_correlations), na.rm = TRUE)
  #0.2018005
}