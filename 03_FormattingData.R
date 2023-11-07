######################Formating data#########################
#######################Date: 2023.10.10######################
######################Author: Yidie##########################

#### Organize exposure_outcome.matchV1.tsv data format into MR Analysis format 
####Add proxy related columns
library(TwoSampleMR)

##Example exposure and outcomes
exposures <- c("LST","MVPA","AbPA","MVPA_sen","ST")
outcomes <- c("AIS","LAS","CAS","SVS")
for (exposure in exposures) {
  for (outcome in outcomes) {
    # read data
    exposure_file <- paste0("./clumped_plink/", exposure, ".clumped.tsv")
    outcome_file <- paste0("./match_V1/", exposure,"_", outcome, ".matchV1.tsv")
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
    
    # Output
    outdat_file <- paste0("./match_V2/", exposure,"_", outcome, ".matchV2.tsv")
    write.table(outcome_dat, file = outdat_file, sep = "\t", quote = FALSE, row.names = FALSE)
  }
}

