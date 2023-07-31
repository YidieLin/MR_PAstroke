###########################创建用于UVMR分析的数据集########################
################Date: 2023.0606###################
###############Author: Yidie ######################
library(TwoSampleMR)
library(plyr)
library(MRPRESSO)
library(mr.raps)

#暴露列表和结局列表
exposures <- c("LST", "MVPA", "AbPA", "BMI", "BFP", "WHR", "SMOKING", "SBP", "LDL_C", "HDL_C","TG")
outcomes <- c("AIS", "CAS", "LAS", "SVS")

for (exposure1 in exposures) {
  for (outcome1 in outcomes) {
    
    print(paste0("BEGIN: 准备数据 ", exposure1,"_", outcome1))
    #prepare data
    #step 1: prepare the exposure data
    exposure <- read.table(paste0("/media/yidie/Linuxdata/clumped/",exposure1,".clumped.tsv"), 
                           header = T, sep = "\t", stringsAsFactors = FALSE)
    
    #orienta to the positive direction
    exposure$EA_c <- exposure$EA; exposure$NEA_c <- exposure$NEA
    exposure[ which(exposure$BETA < 0), "EA"] <-   exposure[ which(exposure$BETA < 0), "NEA_c"]      #effect allele
    exposure[ which(exposure$BETA < 0), "NEA"] <-   exposure[ which(exposure$BETA < 0), "EA_c"]      #non effect allele
    exposure[ which(exposure$BETA < 0), "EAF"] <-  1 - exposure[ which(exposure$BETA < 0), "EAF"]      #effect allele frequency
    exposure$BETA <- abs(exposure$BETA)  
    exposure <- exposure[, c(1:11)] # SNP CHR POS EA NEA EAF BETA SE P. need to check the EAF colunm before
    ##
    exposure_dat <- format_data(exposure, type = "exposure", header = TRUE, phenotype_col = "PHENOTYPE",
                                snp_col = "SNP", beta_col = "BETA", se_col = "SE", eaf_col = "EAF", effect_allele_col = "EA",
                                other_allele_col = "NEA", pval_col = "P", samplesize_col = "N")
    exposure_dat_tmp <- exposure_dat
    print(paste("step 1: 准备暴露数据, 输入", nrow(exposure), "SNPs", "。暴露数据集的维度：", dim(exposure)[1], dim(exposure)[2]))
    ##
    
    #step 2: prepare the outcome data by extract the summary statistics using above IVs information
    outcome_dat <- read.table(paste0("/media/yidie/Workspace/analysis0607_stroke/data/matched/", exposure1, "_", outcome1, ".clumped2.tsv"), 
                              header = T, sep = "\t", stringsAsFactors = FALSE)
    
    print(paste("step 2: 根据暴露IVs的信息，从结局数据集中获取SNPs: ", nrow(outcome_dat), "个。 结局数据集的维度：", dim(outcome_dat)[1], dim(outcome_dat)[2]))
    ######################
    #step 3: remove potential pleiotropic SNPs which is genome-wide significant associated (<5E-8) with other risk factors
    # "AIS", "CAS", "LAS","SVS","BMI", "BFP", "WHR"
    list_func <- function(var, dat){
      if (var %in% c("LST", "MVPA", "AbPA", "BMI", "BFP", "WHR", "SMOKING", "SBP", "LDL_C", "HDL_C","TG"))  dat1 <- dat[(dat$Exposure %in% c("AIS", "LAS", "CAS", "SVS")), ]
      
      return(dat1)
    }
    
    dat_pleio <- read.table(paste0("/media/yidie/Workspace/analysis0607_stroke/results/pleio_signSNP_list"), header = T)
    
    print(paste("the dim of pleiotropy snps is", dim(dat_pleio)[1], dim(dat_pleio)[2]))
    
    dat_pleio1 <- list_func(var = "LST", dat = dat_pleio)
    
    exposure_dat_tmp1 <- exposure_dat_tmp[!(exposure_dat_tmp$SNP %in% dat_pleio1$SNP), ]
    exposure_dat_pleio <- exposure_dat_tmp[(exposure_dat_tmp$SNP %in% dat_pleio1$SNP), ]
    write.table(exposure_dat_pleio, file = paste0("/media/yidie/Workspace/analysis0607_stroke/data/pleioSNPs/",exposure1,"_",outcome1,".pleiosnp"), sep = "\t", quote = FALSE, row.names = FALSE)
    
    print(paste("step 3: remove potential pleiotropic SNPs which is genome-wide significant associated (<5E-8) with other risk factors in current study using the meta's dataset of risk factors, leave", nrow(exposure_dat_tmp1), "SNPs", ", the dim of exposure_dat_tmp1 dataset is", dim(exposure_dat_tmp1)[1], dim(exposure_dat_tmp1)[2]))
    #############################
    
    #step 4: harmonise the exposure and outcome data
    e_o_dat <- harmonise_data(exposure_dat = exposure_dat_tmp1, outcome_dat = outcome_dat, action = 1)
    e_o_dat_clean <- e_o_dat[e_o_dat$mr_keep == "TRUE", ]
    print(paste("step 4: 统一暴露与结局数据, 保留：", nrow(e_o_dat_clean), "SNPs", "。e_o_dat_clean 数据集的维度：", dim(e_o_dat_clean)[1], dim(e_o_dat_clean)[2]))
    ##
    
    #step5-1: 计算Variance explained by IV and IV power，保存剔除outliers前的分析数据集
    e_o_dat_clean$Rsquare <- (2*e_o_dat_clean$eaf.exposure*(1-e_o_dat_clean$eaf.exposure)*(e_o_dat_clean$beta.exposure)^2)
    e_o_dat_clean$Fstat <- (e_o_dat_clean$Rsquare* (e_o_dat_clean$samplesize.exposure-2)/(1-e_o_dat_clean$Rsquare))
    
    #输出Harmonised数据
    outdat_file <- paste0("/media/yidie/Workspace/analysis0607_stroke/data/analysisData/", exposure1, "_", outcome1, ".csv")
    write.table(e_o_dat_clean, file = outdat_file, sep = ",", quote = FALSE, row.names = FALSE)
    ##
    print("step 4: 经过工具变量筛选流程, e_o_dat_clean数据集用于进行后续MR因果主分析")
    
    print("-------------------------------------------------------------------------")
    
  }
  
}



############################################
#########################################
#暴露列表和结局列表
exposures <- c("LST", "MVPA", "AbPA", "BMI", "BFP", "WHR", "SMOKING", "SBP", "LDL_C", "HDL_C","TG")
outcomes <- c("AIS", "CAS", "LAS", "SVS")
#step 5: 采用MR PRESSO test和leave-one-out分析Outliers并保存Outlier数据
for (exposure1 in exposures) {
  for (outcome1 in outcomes) {
    #step 1: prepare the exposure data
    e_o_dat_clean <- read.table(paste0("/media/yidie/Workspace/analysis0607_stroke/data/analysisData/", exposure1, "_", outcome1, ".csv"), 
                                header = T, sep = ",", stringsAsFactors = FALSE)
    print(paste0("数据导入完毕：", exposure1,"_", outcome1))
    
    ##(1)采用PRESSO检测
    if (nrow(e_o_dat_clean) > 5){
      
      #initial
      outlier_index <- "NA"
      dat_presso <- e_o_dat_clean
      
      #while loop
      while ( any(!is.na(outlier_index)) ){
        
        mrpresso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", 
                              OUTLIERtest = TRUE, DISTORTIONtest = TRUE, 
                              data = dat_presso, 
                              NbDistribution = 5000, 
                              SignifThreshold = 0.05, 
                              seed = 20230405)
        
        if (is.null(mrpresso$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`)) {
          outlier_index <- NA
        } else if ("No significant outliers" %in% mrpresso$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`) {
          outlier_index <- NA
        } else {
          outlier_index <- row.names(mrpresso$`MR-PRESSO results`$`Outlier Test`[mrpresso$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`, ])
        }
        
        dat_presso <- dat_presso[!(rownames(dat_presso) %in% c(outlier_index)), ]  
        
      }
      #mark outlier SNP
      e_o_dat_clean$outlier <- "no";
      e_o_dat_clean$outlier[!(rownames(e_o_dat_clean) %in% rownames(dat_presso ))] <- "yes"                               
      
      print(paste("step 5: MR PRESSO检验outliers，发现", nrow(e_o_dat_clean) - nrow(dat_presso), "outlier SNPs，并剩余 ", nrow(dat_presso), "SNPs。dat_presso 的维度：", dim(dat_presso)[1], dim(dat_presso)[2]))
      print(paste0("检测出的Outliers有: ",e_o_dat_clean$SNP[e_o_dat_clean$outlier == "yes"]))
    } else {
      e_o_dat_clean$outlier <- "no";
      dat_presso  <- e_o_dat_clean
      
      print(paste("step 5: MR PRESSO检验outliers，发现0 outlier SNPs，保留", nrow(dat_presso), "SNPs。e_o_dat_clean1的维度：", dim(dat_presso)[1], dim(dat_presso)[2]))
    }
    
    ######(2)采用leave-one-out检测outliers
    res_loo <- mr_leaveoneout(e_o_dat_clean)
    # 获取SNP值为"All"的行的beta值
    all_values <- res_loo[res_loo[, 6] == "All", 7]
    
    # 获取Outliers
    if (any(all_values > 0)) {
      outlier_loo <- res_loo[res_loo[, 7] < 0, 6]
    } else {
      outlier_loo <- res_loo[res_loo[, 7] > 0, 6]
    }
    
    print(paste("step 5: Leave-one-out 发现", nrow(outlier_loo), "个outlier SNP"))
    
    
    ##创建最终的outliers数据集
    outlier_presso <- e_o_dat_clean[e_o_dat_clean$outlier == "yes", 1]
    data_outlier_1 <- data.frame(outlier_presso)
    data_outlier_2 <- data.frame(outlier_loo)
    
    #输出Outliers数据
    outdat_file <- paste0("/media/yidie/Workspace/analysis0607_stroke/data/outliers/", exposure1, "_", outcome1, ".outlierpresso")
    write.table(data_outlier_1, file = outdat_file, sep = ",", quote = FALSE, row.names = FALSE)
    outdat_file2 <- paste0("/media/yidie/Workspace/analysis0607_stroke/data/outliers/", exposure1, "_", outcome1, ".outlierloo")
    write.table(data_outlier_2, file = outdat_file2, sep = ",", quote = FALSE, row.names = FALSE)
    ##
    print(paste("Step5: MRPRESSO已成功检测outlier", data_outlier_1$outlier_presso))
    print(paste("Step5: Leave-one-out 已成功检测outlier", data_outlier_2$outlier_loo))
    print(paste("-------------------------------END：分割线--------------------------------------"))
    
    ##
    gc()
  }
  
}
