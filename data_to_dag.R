#exposure:antidepressants, outcome:white blood cell count (WBC)
#med_tvc denotes whether they were taking antidepressants

df = read.table("~/server/sealockj/projects/mdd_wbc_long/problem_list_meds/time_censor/time_censor_meds_and_labs/pgs_moderation/datasets/mdd_pgs_moderation_between_wbc_and_TCA_365_days_no_inpt_dataset_2021-06-23.txt", 
                header = TRUE, stringsAsFactors = FALSE)
