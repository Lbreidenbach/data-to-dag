#exposure:antidepressants, outcome:white blood cell count (WBC)
#med_tvc denotes whether they were taking antidepressants at time of measurement
library(bnlearn)


wbc_df = read.table("~/server/sealockj/projects/mdd_wbc_long/problem_list_meds/time_censor/time_censor_meds_and_labs/pgs_moderation/datasets/mdd_pgs_moderation_between_wbc_and_TCA_365_days_no_inpt_dataset_2021-06-23.txt", 
                header = TRUE)



reclassify = function(df){
  int_true = sapply(df, is.integer)
  logi_true = sapply(df, is.logical)
  int_ind = as.numeric(unname(which(int_true, useNames = FALSE)))
  logi_ind = as.numeric(unname(which(logi_true, useNames = FALSE)))
  df[int_ind] = lapply(df[int_ind], as.factor)
  df[logi_ind] = lapply(df[logi_ind], as.factor)
  return(df)
}

wbc_df = reclassify(wbc_df)
wbc_df_sub = wbc_df[c(-1, -2, -6:-18, -22, -26, -27, -30, -31)] #cutting out all numeric data for now. figure out how to incorporate numerics later

dag = mmpc(wbc_df_sub, undirected = FALSE)
plot(dag)

#fig

med_1_df = df[df$med_tvc == 1, ]
med_0_df = df[df$med_tvc == 0, ]
