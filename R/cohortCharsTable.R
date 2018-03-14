###################################################################
# Patient characteristics table
###################################################################

cohortchars = function(df, colvar){
   #--------------------------------------------------------------
   # FUNCTION: Creates a Table 1 of cohort characteristics, using
   #           the package qwraps2.  Note: Need to manually update
   #           the summaries object for each new data frame.
   #--------------------------------------------------------------
   # INPUTS:   df = the data frame
   #           colvar = the column variable
   #--------------------------------------------------------------
   # OUTPUT:   a data frame containing a cohort chars table
   #--------------------------------------------------------------
   
   require(rlist)
   require(tidyverse)
   require(qwraps2)
   options(qwraps2_markup = 'markdown')
   
   # Set up summaries info.
   summaries <-
      list("Gestational age (w):" = 
              list("28-31, n (%)"           = ~ qwraps2::n_perc0(gest_age <= 31),
                   "32-36, n (%)"           = ~ qwraps2::n_perc0(gest_age %in% 32:36),
                   "37-38, n (%)"           = ~ qwraps2::n_perc0(gest_age %in% 37:38),
                   "39-40, n (%)"           = ~ qwraps2::n_perc0(gest_age %in% 39:40),
                   "41-43, n (%)"           = ~ qwraps2::n_perc0(gest_age >= 41)
              ),
           "Infant sex:" =
              list("male, n (%)"            = ~ qwraps2::n_perc0(sex=='male'),
                   "female, n (%)"          = ~ qwraps2::n_perc0(sex=='female')
              ),
           "Maternal age (y):" = 
              list("<20, n (%)"             = ~ qwraps2::n_perc0(mat_age < 20),
                   "20-29, n (%)"           = ~ qwraps2::n_perc0(mat_age %in% 20:29),
                   "30-39, n (%)"           = ~ qwraps2::n_perc0(mat_age %in% 30:39),
                   "40-49, n (%)"           = ~ qwraps2::n_perc0(mat_age %in% 40:49),
                   ">=50, n (%)"            = ~ qwraps2::n_perc0(mat_age >= 50)
              ),
           "Maternal job type:" = 
              list("low, n (%)"             = ~ qwraps2::n_perc0(mat_jobtype == 'low'),
                   "medium, n (%)"          = ~ qwraps2::n_perc0(mat_jobtype == 'med'),
                   "high, n (%)"            = ~ qwraps2::n_perc0(mat_jobtype == 'high')
              ),
           "Stillborn:" =
              list("No, n (%)"              = ~ qwraps2::n_perc0(sb==0),
                   "Yes, n (%)"             = ~ qwraps2::n_perc0(sb==1)
              ),
           "Parity:" =
              list("Primiparous, n (%)"     = ~ qwraps2::n_perc0(primip==1),
                   "Multiparous, n (%)"     = ~ qwraps2::n_perc0(primip==0)
              ),
           "HIV:" =
              list("No, n (%)"              = ~ qwraps2::n_perc0(HIV==1),
                   "Yes, n (%)"             = ~ qwraps2::n_perc0(HIV==0)
              ),
           "Blood pressure:" = 
              list('Systolic, mean (sd)'    = ~ qwraps2::mean_sd(bpsys_high, denote_sd = "paren"),
                   'Diastolic, mean (sd)'   = ~ qwraps2::mean_sd(bpdias_high, denote_sd = "paren")
              ), 
           "Protein level:" = 
              list('0, n (%)'               = ~ qwraps2::n_perc0(protein_level==0),
                   '1, n (%)'               = ~ qwraps2::n_perc0(protein_level==1),
                   '2, n (%)'               = ~ qwraps2::n_perc0(protein_level==2),
                   '3, n (%)'               = ~ qwraps2::n_perc0(protein_level==3),
                   '4, n (%)'               = ~ qwraps2::n_perc0(protein_level==4)
              )
      )
   
   # Create table.
   table = cbind(
      summary_table(df, summaries), # For entire df
      summary_table(dplyr::group_by_(df, colvar), summaries)) # For newpet columns
   
   # Demographic names and subgroup names.
   names = sort(c(list.names(summaries),names(list.flatten(summaries))))
   names = gsub(".*\\:\\.","",names)
   
   # Indices for which rows contain data (using :).
   idx = (1:length(names))[-grep(":",names)]
   idx_grps = grep(":",names)
   
   # Convert table to data frame.
   table1 = data.frame('Demographic' = names)
   
   for(i in 1:ncol(table)){
      table1[idx,i+1] = table[,i]
      table1[idx_grps,i+1] = ""
   }
   
   # Set column names.
   colnames(table1)[-1] = colnames(table)
   colnames(table1)[2] = gsub('df','cohort',colnames(table1)[2])
   return(table1)
}


