###########################多变量MR分析########################
################Date: 2023.0608###################
###############Author: Yidie ######################


#######创建PA、肥胖相关形状与卒中的多变量MR分析函数#######
robMV1 <- function(exposure1, exposure2, outcome){
  library(TwoSampleMR)
  library(MVMR)
  library(robustMVMR)
  library(ggplot2)
  library(robustbase)
  library(stats)
  ## data preparation
  # step1: prepare the initial IVs list
  e1_uvdat <- read.delim(paste0("/home/yidie/work_dir/data/summary_data/PA/",exposure1,".clumped.tsv"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  e2_uvdat <- read.delim(paste0("/home/yidie/work_dir/data/summary_data/",exposure2,"/",exposure2,".clumped.tsv"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  e_uvdat <- rbind(e1_uvdat,e2_uvdat)  
  
  # remove duplicated SNPs
  e_uvdat <- e_uvdat[!duplicated(e_uvdat$SNP), ]
  print(paste("step 1: 获取暴露数据(",exposure1,", ", exposure2, ")的初始工具变量", " 总共 ", nrow(e_uvdat), " SNPs", ", e_uvdat dataset的维度为", dim(e_uvdat)[1], dim(e_uvdat)[2]))  
  
  #step 2: select independent (no LD) variants
  ld_list <- subset(e_uvdat)
  ld_list[which(ld_list$SNP %in% e1_uvdat$SNP), ]$P <- 0  
  names(ld_list)[9] <- "pval.exposure"  
  
  e_uvdat_list <- clump_data(ld_list, clump_kb = 1000, clump_r2 = 0.01, clump_p1 = 1, 
                             clump_p2 = 1, pop = "EUR")
  #the full snps list for subsequently MV MR analysis
  e_uvdat1 <- subset(e_uvdat, SNP %in% e_uvdat_list$SNP)
  #leave SNPs
  print(paste("step 2: 进行LD clump后，独立的variants用于后续分析，剩余：", length(unique(e_uvdat1$SNP)), "SNPs"))
  
  #step 3: remove the pleiotropic snps
  list_func <- function(var, dat){
    if (var %in% c("LST", "AbPA", "MVPA"))  
      dat1 <- dat[(dat$Exposure %in% c("AIS", "LAS", "CAS", "SVS")), ]
    return(dat1)
  }
  
  dat_pleio <- read.delim(file = "/home/yidie/work_dir/analysis0607_nopleio/results/pleio_signSNP_list", sep = "\t" , header = T)
  print(paste("the dim of pleiotropy snps is", dim(dat_pleio)[1], dim(dat_pleio)[2]))
  dat_pleio1 <- list_func(var = exposure1, dat = dat_pleio)
  table(dat_pleio1$Exposure)
  
  e_uvdat2 <- e_uvdat1[!(e_uvdat1$SNP %in% dat_pleio1$SNP), ]
  print(paste("step 3: 去除水平多效性SNPS：", length(unique(e_uvdat1$SNP))-length(unique(e_uvdat2$SNP)), "个。","剩余：",length(unique(e_uvdat2$SNP))," SNPs"))
  ##
  #step4: select variants existed in the exposures' GWAS dataset
  e1_mvdat <- read_outcome_data(snps = unique(e_uvdat2$SNP), 
                                filename = paste0("/media/yidie/Linuxdata/GWAS/EUR/",exposure1,".tsv"),
                                sep = "\t", snp_col = "SNP", beta_col = "BETA",
                                se_col = "SE", effect_allele_col = "EA",
                                other_allele_col = "NEA", eaf_col = "EAF",
                                pval_col = "P", phenotype_col = "PHENOTYPE")
  names(e1_mvdat) <- gsub("outcome", "exposure", names(e1_mvdat))
  e2_mvdat <- read_outcome_data(snps = unique(e_uvdat2$SNP), 
                                filename = paste0("/media/yidie/Linuxdata/GWAS/EUR/",exposure2,".tsv"),
                                sep = "\t", snp_col = "SNP", beta_col = "BETA",
                                se_col = "SE", effect_allele_col = "EA",
                                other_allele_col = "NEA", eaf_col = "EAF",
                                pval_col = "P", phenotype_col = "PHENOTYPE")

  
  #orienta to the positive direction for the primary exposure: as suggested by Andrew J. Grant
  e1_mvdat$effect_allele.exposure_c <- e1_mvdat$effect_allele.exposure
  e1_mvdat$other_allele.exposure_c <- e1_mvdat$other_allele.exposure
  e1_mvdat[which(e1_mvdat$beta.exposure < 0), "effect_allele.exposure"] <- e1_mvdat[which(e1_mvdat$beta.exposure < 0), "other_allele.exposure_c"]
  e1_mvdat[which(e1_mvdat$beta.exposure < 0), "other_allele.exposure"] <- e1_mvdat[which(e1_mvdat$beta.exposure < 0), "effect_allele.exposure_c"]
  e1_mvdat[which(e1_mvdat$beta.exposure < 0), "eaf.exposure"] <- 1-e1_mvdat[which(e1_mvdat$beta.exposure < 0), "eaf.exposure"]
  e1_mvdat$beta.exposure <- abs(e1_mvdat$beta.exposure)
  e1_mvdat <- e1_mvdat[, c(1:12)]
  
  #leave SNPs
  print(paste("step 4: 从主暴露（",exposure1,"）GWAS数据集获取合并的SNPs信息：", nrow(e1_mvdat), "SNPs", ", e1_mvdat数据的维度：", dim(e1_mvdat)[1], dim(e1_mvdat)[2]))
  
  #intersection: keep the shared variants among these exposures
  inter_list <- Reduce(intersect, list(e1_mvdat$SNP, e2_mvdat$SNP))
  ## 
  print(paste("step 4: 再从协变量（", exposure2,"）", "的GWAS数据中获取SNPs信息，总共得到：", length(inter_list), "SNPs"))
  ##
  
  #step5: harmonise snps to be all on the same strand
  e2_e4_mvdat <- rbind(subset(e2_mvdat, SNP %in% inter_list))
  e_mvdat_h <- harmonise_data(e1_mvdat, e2_e4_mvdat, action = 1) 
  
  print(paste("step 5: 统一暴露 (",exposure1,", ",exposure2,") 的SNPs正负链, 剩余：", length(unique(e_mvdat_h$SNP)), "SNPs。e_mvdat_h数据集的维度：", dim(e_mvdat_h)[1], dim(e_mvdat_h)[2]))
  ##
  #step 6: format the exposures' dataset
  e_mvdat_h1 <- subset(e_mvdat_h, id.outcome == id.outcome[1], select = c(SNP, exposure, id.exposure, effect_allele.exposure, other_allele.exposure, eaf.exposure, beta.exposure, se.exposure, pval.exposure))
  e_mvdat_h2 <- subset(e_mvdat_h, select = c(SNP, outcome, id.outcome, effect_allele.outcome, other_allele.outcome, eaf.outcome, beta.outcome, se.outcome, pval.outcome))
  names(e_mvdat_h2) <- gsub("outcome", "exposure", names(e_mvdat_h2))
  e_mvdat_h_qc <- rbind(e_mvdat_h1, e_mvdat_h2)
  
  e_mvdat_h_qc[e_mvdat_h_qc$exposure == exposure1, ]$id.exposure <- exposure1
  e_mvdat_h_qc[e_mvdat_h_qc$exposure == exposure2, ]$id.exposure <- exposure2
  
  e_mvdat_h_qc$id.exposure <- factor(e_mvdat_h_qc$id.exposure, levels = c(exposure1, exposure2), 
                                     labels = c(exposure1, exposure2))
  
  print(paste("step 6: 整理并合并暴露数据集 (",exposure1, " 和 ",exposure2, ") 的格式, 保留：", length(unique(e_mvdat_h_qc$SNP)), "SNPs。 e_mvdat_h_qc数据集的维度：", dim(e_mvdat_h_qc)[1], dim(e_mvdat_h_qc)[2]))
##
  
  #step 7: select the variants existed in the outcome' GWAS dataset 
  outcome_mvdat <- read_outcome_data(snps = unique(e_mvdat_h_qc$SNP), 
                                     filename = paste0("/home/yidie/work_dir/analysis0607_nopleio/data/MVMR_matched/",exposure1,"_",exposure2,"_",outcome,".clumped.tsv"), 
                                     sep = "\t", snp_col = "SNP", beta_col = "beta.outcome", se_col = "se.outcome", effect_allele_col = "effect_allele.outcome", 
                                     other_allele_col = "other_allele.outcome", eaf_col = "eaf.outcome", pval_col = "pval.outcome", phenotype_col = "outcome") 
  #leave SNPs
  e_mvdat_h_qc_2 <- e_mvdat_h_qc[(e_mvdat_h_qc$SNP %in% outcome_mvdat$SNP), ]
 
  print(paste("step 7: 获取结局(",outcome,")GWAS数据集中存在的SNPs，保留：", length(unique(e_mvdat_h_qc_2$SNP)), "SNPs"))
  ##  
  
  #step8: harmonise the exposure and outcome data
  mvdat <- mv_harmonise_data(exposure_dat = e_mvdat_h_qc_2, outcome_dat = outcome_mvdat, harmonise_strictness = 1)
  #leave SNPs
  print(paste("step 8: 协调暴露与结局数据后，进行后续MVMR分析，MVMR数据包含", dim(mvdat$exposure_beta)[1], "SNPs"))
  ##  
  ##
  
  #MVMR analysis
  #step9: multivariable MR analysis
  #based on robustMVMR
  #format data
  mvdat2 <- format_mvmr(BXGs = mvdat$exposure_beta, BYG = mvdat$outcome_beta,
                        seBXGs = mvdat$exposure_se, seBYG = mvdat$outcome_se)
  mv_dat <- mvdat2
  mv_dat$SNP <- row.names(mvdat2)
  
  write.table(mv_dat, file = paste0("/home/yidie/work_dir/analysis0607_nopleio/data/analysisData_MVMR1/",exposure1,"_",exposure2,"_",outcome,".tsv"), 
              sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  ##  
  rob_dat <- mvdat2
  rob_dat <- matrix(unlist(mvdat2), nrow(mvdat2))
  dimnames(rob_dat) <- list(c(1:nrow(rob_dat)), c("SNP", "beta.Y", "sebeta.Y", "beta.e1", "beta.e2",
                                                  "sebeta.e1", "sebeta.e2"))
  pvalue.Y <- 5e-8
  pvalue.e1 <- 5e-8
  pvalue.e2 <- 5e-8

  rob_dat1 <- cbind(rob_dat, pvalue.Y, pvalue.e1, pvalue.e2)
  betaGY <- rob_dat1[, "beta.Y"]
  sebetaGY <- rob_dat1[, "sebeta.Y"]; pvalbetaGY <- rob_dat1[, "pvalue.Y"]
  betaGX <- rob_dat1[, c("beta.e1","beta.e2")]
  sebetaGX <- rob_dat1[, c("sebeta.e1","sebeta.e2")]
  pvalbetaGX <- rob_dat1[, c("pvalue.e1","pvalue.e2")]
  
  #estimate effect size
  robMVMR <- robustMVMR(betaGY = betaGY, sebetaGY = sebetaGY, pvalbetaGY = pvalbetaGY,
                        betaGX = betaGX, sebetaGX = sebetaGX, pvalbetaGX = pvalbetaGX,
                        pval_threshold = 1, plot = FALSE)
  
  #output all result
  #robustMVMR
  robmv_res_out1 <- subset(as.data.frame(robMVMR$mvMRResult_heter_robust), select = c(Index, betaXY, sebetaXY, PvalbetaXY, cond_F_stat))
  robmv_res_out1$betaXY <- formatC(robmv_res_out1$betaXY, format = 'f', digits = 4)
  robmv_res_out1$sebetaXY <- formatC(robmv_res_out1$sebetaXY, format = 'f', digits = 4)
  robmv_res_out1$PvalbetaXY <- formatC(robmv_res_out1$PvalbetaXY, format = 'e', digits = 2)
  robmv_res_out1$cond_F_stat <- formatC(robmv_res_out1$cond_F_stat, format = 'f', digits = 4)
  names(robmv_res_out1)[2:5] <- c("Beta", "SE", "P", "Fstatistic")
  robmv_res_out1$Exposure <- c(exposure1, exposure2)
  robmv_res_out1$Method <- "Based on Grant (MVMR-Robust)"
  robmv_res_out1 <- subset(robmv_res_out1, select = c(Method, Index, Exposure, Beta, SE, P, Fstatistic))      
  write.table(robmv_res_out1, file = paste0("/home/yidie/work_dir/analysis0607_nopleio/results/MVMR1/",exposure1,exposure2,"_",outcome,"rMVMR.csv"), 
              append = TRUE, row.names = F, col.names = T, quote = F, sep = ",")
  #
  ##
  print(paste("step 9: 完成 ",exposure1, ", ", exposure2, "对",outcome,"的MVMR分析，并输出结果。"))
  ## 

}

#定义需要处理的exposure1、exposure2和outcome列表
 exposure1_list <- c("LST", "MVPA", "AbPA")
 exposure2_list <- c("BMI", "BFP", "WHR")
 outcome_list <- c("AIS", "LAS", "CAS","SVS")

 
robMV1(exposure1 = "LST", exposure2 = "BMI", outcome = "AIS")
 
 
#使用for循环批量运行robMR1函数
for (i in 1:length(exposure1_list)) {
  for (j in 1:length(exposure2_list)) {
    for (k in 1:length(outcome_list)) {
      # call the robMV1 function with the corresponding exposures and outcome
      robMV1(exposure1 = exposure1_list[i], 
             exposure2 = exposure2_list[j], 
             outcome = outcome_list[k])
    }
  }
}
