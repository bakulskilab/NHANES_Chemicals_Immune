ensure_pop_detect_freq_same_regression <- function(df_inclusion_criteria
                                                   , df_regression)
{
  df_nobs_merged <- df_inclusion_criteria %>%
    left_join(.
              , df_regression %>%
                rename(chemical_codename_use = chemical_codename) %>%
                filter(term == "chem_log_measurement") %>%
                select(chemical_codename_use
                       , nobs))
  print(colnames(df_nobs_merged))
  
  subset_nobs <- df_nobs_merged %>%
    select("chemical_codename_use"      
           , "chemical_name"
           , "total_number_people_unweighted"
           , "nobs") %>%
    mutate(diff = total_number_people_unweighted - nobs)
  View(subset_nobs)
  
  problem_chems <- subset_nobs %>%
    filter(diff != 0) %>%
    pull(chemical_codename_use)
  # print(problem_chems)
  
  return(problem_chems)
}