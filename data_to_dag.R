#exposure:antidepressants, outcome:white blood cell count (WBC)
#med_tvc denotes whether they were taking antidepressants at time of measurement
library(bnlearn)


df = read.table("~/server/sealockj/projects/mdd_wbc_long/problem_list_meds/time_censor/time_censor_meds_and_labs/pgs_moderation/datasets/mdd_pgs_moderation_between_wbc_and_TCA_365_days_no_inpt_dataset_2021-06-23.txt", 
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

df = reclassify(df)
num_true = sapply(df, is.numeric)
num_ind = as.numeric(unname(which(num_true, useNames = FALSE)))
#df[num_ind]

#applies all discretize methods to numeric data
list_df = lapply(X = c("interval", "quantile", "hartemink"), 
                 FUN = function(method) discretize(
                   data = df[num_ind],
                   method = method,
                   breaks = 4,
                   ordered = TRUE
                 ))

names(list_df) = c("interval", "quantile", "hartemink")

#wbc_df_sub = wbc_df[c(-1, -2, -6:-18, -22, -26, -27, -30, -31)] #cutting out all numeric data for now. figure out how to incorporate numerics later
#v_algorithms = c("pc.stable", "gs", "iamb", "fast.iamb", "inter.iamb", "iamb.fdr", "mmpc", 
#                 "si.hilton.pc", "hpc", "hc", "tabu", "rsmax2", "mmhc", "h2pc", "aracne", "chow.liu")
#slower tests, iamb (gs breaks R), hpc and h2pc takes a few seconds, hc doesn't work
v_algorithms = c("fast.iamb", "mmpc", "si.hiton.pc", "pc.stable", "iamb.fdr", "hpc", "hc", "tabu", "rsmax2", "mmhc", "h2pc", "aracne", "chow.liu")

list_bnlearn = list()

# for(j in v_algorithms) for(k in names(list_df))try({
#   list_bnlearn [[j]][[k]] = do.call(
#     what = j,
#     args = list(x = list_df[[k]])
#   )
# })
for(j in v_algorithms) for(k in names(list_df)) try({
  list_bnlearn[[j]][[k]] <- do.call(
    what = j,
    args = list(x = list_df[[k]])
  )
  M_arcs <- arcs(list_bnlearn[[j]][[k]])
  for(l in 1:nrow(M_arcs)){
    list_bnlearn[[j]][[k]] <- set.arc(
      x = list_bnlearn[[j]][[k]],
      from = M_arcs[l,1],
      to = M_arcs[l,2],
      check.cycles = FALSE,
      check.illegal = FALSE
    )
    list_bnlearn[[j]][[k]] <- choose.direction(
      x = list_bnlearn[[j]][[k]],
      arc = M_arcs[l,],
      data = list_df[[k]]
    )
  }
},silent = TRUE)
# df_arcs = arcs(list_bnlearn[[j]][[k]])
# for(l in 1:row(df_arcs)){
#   list_bnlearn[[j]][[k]] = set.arc(
#     
#   )
# }
######Scoring
M_score <- matrix(
  data = NA,
  nrow = length(v_algorithms),
  ncol = length(list_df),
)
rownames(M_score) <- v_algorithms
colnames(M_score) <- names(list_df)

for(j in v_algorithms) for(k in names(list_df)) try({
  M_score[j,k] <- score(
    x = list_bnlearn[[j]][[k]],
    data = list_M[[k]],
    type = "bic"
  )
})
for(j in rownames(M_score)) M_score <- M_score[,order(M_score[j,])]
for(j in colnames(M_score)) M_score <- M_score[order(M_score[,j]),]
M_score
#figure out scoring tommorrow
graphviz.plot(
  list_bnlearn[[nrow(M_score)]][[ncol(M_score)]]
)
# dag = mmpc(wbc_df_sub, undirected = FALSE, debug = TRUE)
# plot(dag)
# 
# #fig
# 
# med_1_df = df[df$med_tvc == 1, ]
# med_0_df = df[df$med_tvc == 0, ]


