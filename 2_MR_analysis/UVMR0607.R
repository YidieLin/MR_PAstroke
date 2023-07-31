###########################UVMR分析########################
################Date: 2023.0606###################
###############Author: Yidie ######################


library(TwoSampleMR)
library(MRPRESSO)
library(mr.raps)

#暴露列表和结局列表
exposures <- c("LST", "MVPA", "AbPA", "BMI", "BFP", "WHR", "SMOKING", "SBP", "LDL_C", "HDL_C","TG")
outcomes <- c("AIS", "CAS", "LAS", "SVS")

####（一）主分析：采用未剔除潜在PRESSO和leave-one-out检验出的Outliers的数据集
for (exposure in exposures) {
  for (outcome in outcomes) {
    # 读取和格式化数据
    mydat <- read.delim(file = paste0("/media/yidie/Workspace/analysis0607_stroke/data/analysisData/", exposure, "_", outcome, ".csv"), header = TRUE, sep = ",")
    
    print(paste0("step1数据读取完毕：", exposure,"_",outcome))
    ####MR分析
    # 进行MR分析(主分析IVW，敏感性分析Egger\Weighted Median\Weighted Mode)
    res <-generate_odds_ratios(mr(mydat, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median")))
    
    print(paste0("step2 mr主分析完毕"))
    # 采用mr.raps包进行基于Robust Adjusted Profile Score的MR分析
    raps <- mr.raps(mydat$beta.exposure, mydat$beta.outcome, mydat$se.exposure, mydat$se.outcome, 
                    over.dispersion = TRUE, 
                    loss.function = "tukey", 
                    diagnosis = FALSE)
    print(paste0("step3 mr-raps分析完毕"))
    ## Mendelian Randomization包
    #进行IVW分析考虑SNP间的相关性，并进行矫正且采用稳健估计
    mydat2 <- dat_to_MRInput(mydat, get_correlations = TRUE)
    IVWcorrel <- MendelianRandomization::mr_ivw(mydat2[[1]], correl = TRUE)
    print(paste0("step4 mr-IVW with correlation分析完毕"))
    
    # 进行Contamixation Mixture method分析
    conmix <- MendelianRandomization::mr_conmix(mydat2[[1]])
    print(paste0("step5 mr-contamixation Mixture method分析完毕"))
    
    ## 其他敏感性分析
    # 异质性检验
    heter <- mr_heterogeneity(mydat)
    # 增加I^2值
    heter$I2 <- paste0(round((heter$Q - heter$Q_df)*100/heter$Q, 2), "%")
    
    # 多效性检验（Egger Intercept）
    pleio <- mr_pleiotropy_test(mydat)
    
    ## 整理结果
    out_res <- subset(res, select = c("outcome", "exposure", "method", "nsnp", "b", 
                                      "se", "pval", "lo_ci", "up_ci", "or", "or_lci95", "or_uci95"))
    out_raps <- data.frame(outcome="outcome", exposure="exposure", method="Robust adjusted profile score", 
                           nsnp=nrow(mydat),
                           b=raps$beta.hat, se=raps$beta.se,pval=raps$beta.p.value, lo_ci=raps$beta.hat-raps$beta.se*1.96,
                           up_ci=raps$beta.hat+raps$beta.se*1.96, or=exp(raps$beta.hat), or_lci95=exp(raps$beta.hat-raps$beta.se*1.96),
                           or_uci95=exp(raps$beta.hat+raps$beta.se*1.96))
    out_IVWcorrel <- data.frame(outcome=IVWcorrel@Outcome, exposure=IVWcorrel@Exposure, method="IVW (incorporate correlation between variants)",
                                nsnp=IVWcorrel@SNPs, b=IVWcorrel@Estimate, se=IVWcorrel@StdError,
                                pval=IVWcorrel@Pvalue, lo_ci=IVWcorrel@CILower, up_ci=IVWcorrel@CIUpper,
                                or=exp(IVWcorrel@Estimate), or_lci95=exp(IVWcorrel@CILower),
                                or_uci95=exp(IVWcorrel@CIUpper))
    out_conmix <- data.frame(outcome=conmix@Outcome, exposure=conmix@Exposure, method="Contamination mixture method",
                             nsnp=conmix@SNPs, b=conmix@Estimate, se="NA",
                             pval=conmix@Pvalue, lo_ci=conmix@CILower, up_ci=conmix@CIUpper,
                             or=exp(conmix@Estimate), or_lci95=exp(conmix@CILower),
                             or_uci95=exp(conmix@CIUpper))
    out_mainres <- rbind(out_res,out_raps,out_IVWcorrel,out_conmix)
    
    out_heter <- subset(heter, select = c("outcome", "exposure", "method", "Q", "Q_df", "Q_pval", "I2"))
    
    out_pleio <- subset(pleio, select = c("outcome", "exposure", "egger_intercept", "se", "pval"))
    
    
    # 输出结果到文件
    output_file_res <- paste0("/media/yidie/Workspace/analysis0607_stroke/results/UVMR/", exposure, "_", outcome, ".res")
    output_file_heter <- paste0("/media/yidie/Workspace/analysis0607_stroke/results/UVMR/", exposure, "_", outcome, ".heter")
    output_file_pleio <- paste0("/media/yidie/Workspace/analysis0607_stroke/results/UVMR/", exposure, "_", outcome, ".pleio")
    
    
    write.table(out_mainres, file = output_file_res, sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(out_heter, file = output_file_heter, sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(out_pleio, file = output_file_pleio, sep = "\t", quote = FALSE, row.names = FALSE)
    
    
    rm("conmix", "heter", "IVWcorrel", "mydat", "mydat2",
       "out_conmix", "out_heter", "out_IVWcorrel", "out_mainres", "out_pleio",
       "out_raps", "out_res",
       "output_file_heter", "output_file_pleio", "output_file_res",
       "pleio", "raps", "res")
    print(paste0("step6: ",exposure,"_", outcome,"结果输出完毕"))
    print(paste0("---------------------------分割线---------------------------"))
  }
}



######(二)剔除outliers并进行重复分析
for (exposure in exposures) {
  for (outcome in outcomes) {
    # 读取和格式化数据
    mydat <- read.delim(file = paste0("/home/yidie/work_dir/analysis0605/data/analysisData/de_outliers/", exposure, "_", outcome, ".csv"), header = TRUE, sep = ",")
    
    ####MR分析
    # 进行MR分析(主分析IVW，敏感性分析Egger\Weighted Median\Weighted Mode)
    res <-generate_odds_ratios(mr(mydat, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median")))
    
    # 采用mr.raps包进行基于Robust Adjusted Profile Score的MR分析
    raps <- mr.raps(mydat$beta.exposure, mydat$beta.outcome, mydat$se.exposure, mydat$se.outcome, 
                    over.dispersion = TRUE, 
                    loss.function = "tukey", 
                    diagnosis = FALSE)
    
    ## Mendelian Randomization包
    #进行IVW分析考虑SNP间的相关性，并进行矫正且采用稳健估计
    mydat2 <- dat_to_MRInput(mydat, get_correlations = TRUE)
    IVWcorrel <- MendelianRandomization::mr_ivw(mydat2[[1]], correl = TRUE)
    
    # 进行Contamixation Mixture method分析
    conmix <- MendelianRandomization::mr_conmix(mydat2[[1]])
    
    ## 其他敏感性分析
    # 异质性检验
    heter <- mr_heterogeneity(mydat)
    # 增加I^2值
    heter$I2 <- paste0(round((heter$Q - heter$Q_df)*100/heter$Q, 2), "%")
    
    # 多效性检验（Egger Intercept）
    pleio <- mr_pleiotropy_test(mydat)
    
    ## 整理结果
    out_res <- subset(res, select = c("outcome", "exposure", "method", "nsnp", "b", 
                                      "se", "pval", "lo_ci", "up_ci", "or", "or_lci95", "or_uci95"))
    out_raps <- data.frame(outcome="outcome", exposure="exposure", method="Robust adjusted profile score", 
                           nsnp=nrow(mydat),
                           b=raps$beta.hat, se=raps$beta.se,pval=raps$beta.p.value, lo_ci=raps$beta.hat-raps$beta.se*1.96,
                           up_ci=raps$beta.hat+raps$beta.se*1.96, or=exp(raps$beta.hat), or_lci95=exp(raps$beta.hat-raps$beta.se*1.96),
                           or_uci95=exp(raps$beta.hat+raps$beta.se*1.96))
    out_IVWcorrel <- data.frame(outcome=IVWcorrel@Outcome, exposure=IVWcorrel@Exposure, method="IVW (incorporate correlation between variants)",
                                nsnp=IVWcorrel@SNPs, b=IVWcorrel@Estimate, se=IVWcorrel@StdError,
                                pval=IVWcorrel@Pvalue, lo_ci=IVWcorrel@CILower, up_ci=IVWcorrel@CIUpper,
                                or=exp(IVWcorrel@Estimate), or_lci95=exp(IVWcorrel@CILower),
                                or_uci95=exp(IVWcorrel@CIUpper))
    out_conmix <- data.frame(outcome=conmix@Outcome, exposure=conmix@Exposure, method="Contamination mixture method",
                             nsnp=conmix@SNPs, b=conmix@Estimate, se="NA",
                             pval=conmix@Pvalue, lo_ci=conmix@CILower, up_ci=conmix@CIUpper,
                             or=exp(conmix@Estimate), or_lci95=exp(conmix@CILower),
                             or_uci95=exp(conmix@CIUpper))
    out_mainres <- rbind(out_res,out_raps,out_IVWcorrel,out_conmix)
    
    out_heter <- subset(heter, select = c("outcome", "exposure", "method", "Q", "Q_df", "Q_pval", "I2"))
    
    out_pleio <- subset(pleio, select = c("outcome", "exposure", "egger_intercept", "se", "pval"))
    
    
    # 输出结果到文件
    output_file_res <- paste0("/media/yidie/Workspace/analysis0607_stroke/results/UVMR/", exposure, "_", outcome, ".outres")
    output_file_heter <- paste0("/media/yidie/Workspace/analysis0607_stroke/results/UVMR/", exposure, "_", outcome, ".outheter")
    output_file_pleio <- paste0("/media/yidie/Workspace/analysis0607_stroke/results/UVMR/", exposure, "_", outcome, ".outpleio")
    
    
    write.table(out_mainres, file = output_file_res, sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(out_heter, file = output_file_heter, sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(out_pleio, file = output_file_pleio, sep = "\t", quote = FALSE, row.names = FALSE)
    
    
    rm("conmix", "heter", "IVWcorrel", "mydat", "mydat2",
       "out_conmix", "out_heter", "out_IVWcorrel", "out_mainres", "out_pleio",
       "out_raps", "out_res", "outcome",
       "output_file_heter", "output_file_pleio", "output_file_res",
       "pleio", "raps", "res")
  }
}

