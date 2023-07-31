library('TwoSampleMR')

#暴露列表和结局列表
exposures <- c("LST", "MVPA", "AbPA", "BMI", "BFP", "WHR", "SMOKING", "SBP", "LDL_C", "HDL_C","TG")
outcomes <- c("AIS", "CAS", "LAS", "SVS")

for (exposure in exposures) {
  for (outcome in outcomes) {
    # 读取和格式化数据
    exposure_file <- paste0("/media/yidie/Linuxdata/clumped/", exposure, ".clumped.tsv")
    outcome_file <- paste0("/media/yidie/Workspace/analysis0607_stroke/data/matched/", exposure,"_", outcome, ".clumped.tsv")
    exposure_dat <- read.delim(file = exposure_file, header = TRUE, sep = "\t")
    outcome_dat <- read.delim(file = outcome_file, header = TRUE, sep = "\t")
    exposure_dat <- format_data(exposure_dat, 
                                type = "exposure", 
                                header = TRUE,
                                phenotype_col = "PHENOTYPE",
                                snp_col = "SNP",
                                beta_col = "BETA",
                                se_col = "SE",
                                eaf_col = "EAF",
                                effect_allele_col = "EA",
                                other_allele_col = "NEA",
                                pval_col = "P",
                                samplesize_col = "N"
    )
    outcome_dat <- format_data(outcome_dat, 
                               type = "outcome", 
                               header = TRUE,
                               snps = exposure_dat$SNP, 
                               phenotype_col = "PHENOTYPE",
                               snp_col = "SNP",
                               beta_col = "BETA",
                               se_col = "SE",
                               eaf_col = "EAF",
                               effect_allele_col = "EA",
                               other_allele_col = "NEA",
                               pval_col = "P",
                               samplesize_col = "N"   
    )
    outcome_dat$id.outcome <- outcome
    outcome_dat$proxy.outcome <- NA
    outcome_dat$target_snp.outcome <- NA
    outcome_dat$proxy_snp.outcome <- NA
    outcome_dat$target_a1.outcome <- NA
    outcome_dat$target_a2.outcome <- NA
    outcome_dat$proxy_a1.outcome <- NA
    outcome_dat$proxy_a2.outcome <- NA
    
    # 输出数据
    outdat_file <- paste0("/media/yidie/Workspace/analysis0607_stroke/data/matched/", exposure,"_", outcome, ".matched.tsv")
    write.table(outcome_dat, file = outdat_file, sep = "\t", quote = FALSE, row.names = FALSE)
  }
}

