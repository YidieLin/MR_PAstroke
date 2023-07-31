###########################创建用于MVMR分析的结果数据集########################
################Date: 2023.0608###################
###############Author: Yidie ######################

####
MV_dat <- function(exposure1, exposure2, outcome){
e1_o_dat <- read.delim(paste0("/home/yidie/work_dir/analysis0607_nopleio/data/matched/",exposure1,"_",outcome,".clumped2.tsv"), 
                       header = TRUE, sep = "\t", stringsAsFactors = FALSE)
e2_o_dat <- read.delim(paste0("/home/yidie/work_dir/analysis0607_nopleio/data/matched/",exposure2,"_",outcome,".clumped2.tsv"), 
                       header = TRUE, sep = "\t", stringsAsFactors = FALSE)
mv_dat <- rbind(e1_o_dat, e2_o_dat)

write.table(mv_dat, file = paste0("/home/yidie/work_dir/analysis0607_nopleio/data/MVMR_matched/",exposure1,"_",exposure2,"_",outcome,".clumped.tsv"), 
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)


}


exposure1_list <- c("LST", "MVPA", "AbPA")
exposure2_list <- c("BMI", "BFP", "WHR")
outcome_list <- c("AIS", "LAS", "CAS","SVS")

#使用for循环批量运行MV_dat函数
for (i in 1:length(exposure1_list)) {
  for (j in 1:length(exposure2_list)) {
    for (k in 1:length(outcome_list)) {
      # call the robMV1 function with the corresponding exposures and outcome
      MV_dat(exposure1 = exposure1_list[i], 
             exposure2 = exposure2_list[j], 
             outcome = outcome_list[k])
    }
  }
}
