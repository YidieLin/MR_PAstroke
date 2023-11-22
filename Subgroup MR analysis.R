##############################################################
###########################Subgroup MR########################
#########################Date: 2023.1021####################
###########################Author: Yidie#####################
######MEGASTROKE
#######AIS==ebi-a-GCST006908
#######CAS==ebi-a-GCST006910
#######SVS==ebi-a-GCST006909
#######LAS==ebi-a-GCST006907

###########################Part1 Create datasets for analysis########################

library(TwoSampleMR)
library(plyr)
library(MRPRESSO)
library(mr.raps)
library(tidyverse)

exposures <- c("LST","MVPA","AbPA","MVPA_2018","ST")
outcomes <- c("ebi-a-GCST006908","ebi-a-GCST006910","ebi-a-GCST006907","ebi-a-GCST006909")

### Prepare data for MR
for (exposure1 in exposures) {
  for (outcome1 in outcomes) {
    
    print(paste0("BEGIN: prepare data ", exposure1,"_", outcome1))
    #prepare data
    #step 1: prepare the exposure data
    exposure <- read.table(paste0("./clumped_plink/",exposure1,".clumped.tsv"), 
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
    print(paste("step 1. ", nrow(exposure), "SNPs", ". Data dimension: ", dim(exposure)[1], dim(exposure)[2]))
    ##
    
    #step 2: prepare the outcome data by extract the summary statistics using above IVs information
    outcome_dat <- extract_outcome_data(
      exposure_dat_tmp$SNP,
      outcome1,
      proxies = TRUE,
      rsq = 0.8,
      align_alleles = 1,
      palindromes = 1,
      maf_threshold = 0.3,
      access_token = ieugwasr::check_access_token(),
      splitsize = 10000,
      proxy_splitsize = 500
    )
    print(paste("step 2. SNPs extracted from outcome datasets: ", nrow(outcome_dat), ". Outcome_dat dimension：", dim(outcome_dat)[1], dim(outcome_dat)[2]))
    ######################
    
    
    #step 3: harmonise the exposure and outcome data
    e_o_dat <- harmonise_data(exposure_dat = exposure_dat_tmp, outcome_dat = outcome_dat, action = 1)
    e_o_dat_clean <- e_o_dat[e_o_dat$mr_keep == "TRUE", ]
    print(paste("step 3. Harmonise the exposure and outcome. Leaving ", nrow(e_o_dat_clean), "SNPs", "e_o_dat_clean dimension：", dim(e_o_dat_clean)[1], dim(e_o_dat_clean)[2]))
    ##
    
    #step 4: Variance explained by IV and IV power
    e_o_dat_clean$Rsquare <- (2*e_o_dat_clean$eaf.exposure*(1-e_o_dat_clean$eaf.exposure)*(e_o_dat_clean$beta.exposure)^2)
    e_o_dat_clean$Fstat <- (e_o_dat_clean$Rsquare* (e_o_dat_clean$samplesize.exposure-2)/(1-e_o_dat_clean$Rsquare))
    
    #Output Harmonised data
    outdat_file <- paste0("./Subgroup_MR/data_beforeoutliers/", exposure1, "_", outcome1, "_V1.csv")
    write.table(e_o_dat_clean, file = outdat_file, sep = ",", quote = FALSE, row.names = FALSE)
    ##
    print("step 4: e_o_dat_clean prepared for two-sample MR analysis")
    
    print("------------------------------------ END --------------------------------------")
    
  }
  
}


## Removing SNPs related to confounders and ambiguous SNPs
exposures <- c("LST","MVPA","AbPA","MVPA_2018","ST")
outcomes <- c("ebi-a-GCST006908","ebi-a-GCST006910","ebi-a-GCST006907","ebi-a-GCST006909")

for (outcome in outcomes) {
  for (exposure in exposures) {
    exposure_file <- paste("./phenoscanner_SNP", "/", exposure, "_pleio.txt", sep = "")
    input_file <- paste("./Subgroup_MR/data_beforeoutliers/", exposure, "_", outcome, "_V1.csv", sep = "")
    output_file <- paste("./Subgroup_MR/data_clean/", exposure, "_", outcome, "_V2.csv", sep = "")
    
    exposure_lines <- readLines(exposure_file)
    input_data <- read.csv(input_file)
    
    matching_rows <- input_data[!input_data$SNP %in% exposure_lines, ]
    cleaned_data <- matching_rows[matching_rows$ambiguous != "TRUE", ]
    
    write.csv(cleaned_data, output_file, row.names = FALSE)
  }
}



############ Part 2.  Subgroup MR analysis########################

data_res <- data.frame()
data_heter <- data.frame()
data_egger <- data.frame()
#exposures and outcomes
exposures <- c("LST","MVPA","AbPA","MVPA_2018","ST")
outcomes <- c("ebi-a-GCST006908","ebi-a-GCST006910","ebi-a-GCST006907","ebi-a-GCST006909")

#### Primary analysis
for (exposure in exposures) {
  for (outcome in outcomes) {
    # read formated dataset
    mydat <- read.delim(file = paste0("./Subgroup_MR/data_clean/", exposure, "_", outcome, "_V2.csv"), header = TRUE, sep = ",")
    
    print(paste0("step1 loading data for ", exposure,"_",outcome))
    #### MR analysis
    # IVW + MR Egger\Weighted Median
    res <-generate_odds_ratios(mr(mydat, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median")))
    
    print(paste0("step2 primary MR analysis: DONE"))
    # Robust Adjusted Profile Score MR
    raps <- mr.raps(mydat$beta.exposure, mydat$beta.outcome, mydat$se.exposure, mydat$se.outcome, 
                    over.dispersion = FALSE, 
                    loss.function = "tukey", 
                    diagnosis = FALSE)
    print(paste0("step3 mr-raps analysis: DONE"))
    ## Mendelian Randomization package
    # Contamixation Mixture method
    mydat2 <- dat_to_MRInput(mydat)
    conmix <- MendelianRandomization::mr_conmix(mydat2[[1]])
    print(paste0("step4 mr-contamixation Mixture method: DONE"))
    
    ## Other sensitivity analysis
    # heterogeneity
    heter <- mr_heterogeneity(mydat)
    # Calculaing I^2
    heter$I2 <- NA
    heter$I2 <- ifelse(heter$Q < heter$Q_df, 0, paste0(abs(round((heter$Q - heter$Q_df) * 100 / heter$Q, 2)), "%"))
    
    
    # MR-Egger Intercept
    pleio <- mr_pleiotropy_test(mydat)
    
    ## Tidy results
    out_res <- subset(res, select = c("outcome", "exposure", "method", "nsnp", "b", 
                                      "se", "pval", "lo_ci", "up_ci", "or", "or_lci95", "or_uci95"))
    out_raps <- data.frame(outcome="outcome", exposure="exposure", method="Robust adjusted profile score", 
                           nsnp=nrow(mydat),
                           b=raps$beta.hat, se=raps$beta.se,pval=raps$beta.p.value, lo_ci=raps$beta.hat-raps$beta.se*1.96,
                           up_ci=raps$beta.hat+raps$beta.se*1.96, or=exp(raps$beta.hat), or_lci95=exp(raps$beta.hat-raps$beta.se*1.96),
                           or_uci95=exp(raps$beta.hat+raps$beta.se*1.96))
    out_conmix <- data.frame(outcome=conmix@Outcome, exposure=conmix@Exposure, method="Contamination mixture method",
                             nsnp=conmix@SNPs, b=conmix@Estimate, se="NA",
                             pval=conmix@Pvalue, lo_ci=conmix@CILower, up_ci=conmix@CIUpper,
                             or=exp(conmix@Estimate), or_lci95=exp(conmix@CILower),
                             or_uci95=exp(conmix@CIUpper))
    out_mainres <- rbind(out_res,out_raps,out_conmix)
    
    out_heter <- subset(heter, select = c("outcome", "exposure", "method", "Q", "Q_df", "Q_pval", "I2"))
    
    out_pleio <- subset(pleio, select = c("outcome", "exposure", "egger_intercept", "se", "pval"))
    
    exp_name <- switch(exposure,
                       "LST" = "Leisure screen time",
                       "AbPA" = "Accelerometer-based physical activity",
                       "MVPA" = "Leisure time moderate-to-vigorous physical activity",
                       "MVPA_2018" = "Overall moderate-to-vigorous physical activity",
                       "ST" = "Overall sedentary duration",
                       "No-defined")
    out_name <- switch(outcome,
                       "ebi-a-GCST006908" = "Total ischemic stroke",
                       "ebi-a-GCST006910" = "Cardioembolic stroke",
                       "ebi-a-GCST006907" = "Large artery stroke",
                       "ebi-a-GCST006909" = "Small vessel stroke",
                       "No-defined" )

    # Tidy results
    out_heter$outcome <- out_name
    out_pleio$outcome <- out_name
    out_heter$exposure <- exp_name
    out_pleio$exposure <- exp_name

    out_mainres %>%
      mutate(
        heter_p = ifelse(method == "Inverse variance weighted", out_heter$Q_pval[out_heter$method == "Inverse variance weighted"], NA_real_),
        eggerint_p = ifelse(method == "Inverse variance weighted", out_pleio$pval, NA_real_),
        se = as.numeric(se),
        `Outcomes` = out_name,
        `Exposures` = exp_name,
        `METHOD` = method,
        `SNPs` = nsnp,
        `Beta` = sprintf("%.3f", b),
        `SE` = if_else(is.na(se), "", sprintf("%.5s", se)),
        `OR (95% Cl)` = sprintf("%.2f(%.2f to %.2f)", or, or_lci95, or_uci95),
        `P` = sprintf("%.3e", pval),
        `Q statistic P-value` = if_else(is.na(heter_p), "", sprintf("%.3e", heter_p)),
        `MR-Egger pleiotropy test P-value` = if_else(is.na(eggerint_p), "", sprintf("%.3e", eggerint_p))
      ) -> out_mainres
    

    # Combined all results
    data_res <- rbind(data_res, out_mainres)
    data_heter <- rbind(data_heter, out_heter)
    data_egger <- rbind(data_egger, out_pleio)
    
    print(paste0("step6. ",exposure,"_", outcome," analysis: DONE"))
    print(paste0("----------------------------------------------------------"))
  }
}

data_res <- data_res %>% distinct(Outcomes, Exposures, METHOD, .keep_all =TRUE)

output_file_res <- paste0("./Subgroup_MR/results_main/Results_PAtoStroke.csv")
output_file_heter <- paste0("./Subgroup_MR/results_main/Heterogeneity.csv")
output_file_pleio <- paste0("./Subgroup_MR/results_main/MRegger_intercept.csv")

write.table(data_res, file = output_file_res, sep = ",", quote = FALSE, row.names = FALSE)
write.table(data_heter, file = output_file_heter, sep = ",", quote = FALSE, row.names = FALSE)
write.table(data_egger, file = output_file_pleio, sep = ",", quote = FALSE, row.names = FALSE)



#####################################################################
##################Part 3. Detecting outliers#######################
data_mrpresso <- data.frame()
data_outliers <- data.frame()
data_distort <- data.frame()

#exposures and outcomes
exposures <- c("LST","MVPA","AbPA","MVPA_2018","ST")
outcomes <- c("ebi-a-GCST006908","ebi-a-GCST006910","ebi-a-GCST006907","ebi-a-GCST006909")

#MR PRESSO test
for (exposure1 in exposures) {
  for (outcome1 in outcomes) {
    #step 1: prepare the exposure data
    e_o_dat_clean <- read.table(paste0("./Subgroup_MR/data_clean/", exposure1, "_", outcome1, "_V2.csv"), 
                                header = T, sep = ",", stringsAsFactors = FALSE)
    print(paste0("Loading datasets：", exposure1,"_", outcome1,". Time：", Sys.time ()))
    
    exp_name <- switch(exposure1,
                       "LST" = "Leisure screen time",
                       "AbPA" = "Accelerometer-based physical activity",
                       "MVPA" = "Leisure time moderate-to-vigorous physical activity",
                       "MVPA_2018" = "Overall moderate-to-vigorous physical activity",
                       "ST" = "Overall sedentary duration",
                       "No-defined")
    out_name <- switch(outcome,
                       "ebi-a-GCST006908" = "Total ischemic stroke",
                       "ebi-a-GCST006910" = "Cardioembolic stroke",
                       "ebi-a-GCST006907" = "Large artery stroke",
                       "ebi-a-GCST006909" = "Small vessel stroke",
                       "No-defined" )
    ##(1)MR PRESSO
    if (nrow(e_o_dat_clean) > 5){
      
      #initial
      outlier_index <- "NA"
      dat_presso <- e_o_dat_clean
      
      #while loop
      mrpresso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", 
                            OUTLIERtest = TRUE, DISTORTIONtest = TRUE, 
                            data = e_o_dat_clean, 
                            NbDistribution = 3000, 
                            SignifThreshold = 0.05, 
                            seed = 20230405)
      
      if (is.null(mrpresso$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`)) {
        outlier_index <- NA
        
      } else if ("No significant outliers" %in% mrpresso$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`) {
        outlier_index <- NA
        
      } else {
        outlier_index <- row.names(mrpresso$`MR-PRESSO results`$`Outlier Test`[mrpresso$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`, ])
        out_distort <- as.data.frame(cbind(mrpresso$`MR-PRESSO results`$`Distortion Test`$`Distortion Coefficient`, mrpresso$`MR-PRESSO results`$`Distortion Test`$Pvalue))
        out_distort$exposue <- exp_name 
        out_distort$outcome <- out_name
        out_distort <- out_distort %>% select(3:4,1:2)
        
        data_distort <- rbind(data_distort, out_distort)
      }
      
      dat_presso <- dat_presso[!(rownames(dat_presso) %in% c(outlier_index)), ]  
      
      #mark outlier SNP
      e_o_dat_clean$outlier <- "no";
      e_o_dat_clean$outlier[!(rownames(e_o_dat_clean) %in% rownames(dat_presso ))] <- "yes"                               
      
      print(paste("MR PRESSO detected: ", nrow(e_o_dat_clean) - nrow(dat_presso), "outlier SNPs. Leaving ", nrow(dat_presso), "SNPs. dat_presso dimension：", dim(dat_presso)[1], dim(dat_presso)[2]))
      print(paste0("Outliers were: ",e_o_dat_clean$SNP[e_o_dat_clean$outlier == "yes"]))
    } else {
      e_o_dat_clean$outlier <- "no";
      dat_presso  <- e_o_dat_clean
      
      print(paste("MR PRESSO detected: 0 outlier SNPs. Leaving ", nrow(dat_presso), "SNPs. e_o_dat_clean dimension：", dim(dat_presso)[1], dim(dat_presso)[2]))
    }
    
    out_mrpressomain <- as.data.frame(mrpresso$`Main MR results`)
    out_global <- rbind(mrpresso$`MR-PRESSO results`$`Global Test`[1],mrpresso$`MR-PRESSO results`$`Global Test`[2])
    out_mrpresso <- cbind(out_mrpressomain, out_global)
    out_mrpresso$exposue <- exp_name 
    out_mrpresso$outcome <- out_name
    out_mrpresso <- out_mrpresso %>% select(8:9,2:7)
    
    if ((nrow(e_o_dat_clean) - nrow(dat_presso)) == 0) {
      outliers <- NA
      out_outliers <- as.data.frame(outliers)
      out_outliers$exposue <- exp_name 
      out_outliers$outcome <- out_name
      out_outliers <- out_outliers %>% select(2:3,1)
    } else {
      outliers <- e_o_dat_clean$SNP[e_o_dat_clean$outlier == "yes"]
      out_outliers <- as.data.frame(outliers)
      out_outliers$exposue <- exp_name 
      out_outliers$outcome <- out_name
      out_outliers <- out_outliers %>% select(2:3,1)
    }
    
    
    #Combined results
    data_mrpresso <- rbind(data_mrpresso, out_mrpresso)
    data_outliers <- rbind(data_outliers, out_outliers)
    
    
    #Output datastes that marked outliers
    outdat_file <- paste0("./SubrgoupMR/data_doutliers/", exposure1, "_", outcome1, "_V3.csv")
    write.table(e_o_dat_clean, file = outdat_file, sep = ",", quote = FALSE, row.names = FALSE)
    ##
    print(paste("END: MRPRESSO detected ", exposure1, "-", outcome1," outliers. ", e_o_dat_clean$SNP[e_o_dat_clean$outlier == "yes"]))
    print(paste("---------------------------------------------------------------------------"))
    print(paste("DONE：", Sys.time ()))
    ##
    gc()
  }
  
}

data_mrpresso <- cbind(data_mrpresso[1:7], unlist(data_mrpresso$RSSobs))
outdat_file_mrpresso <- paste0("./Subgroup_MR/sensitivity/MRPRESSO.csv")
outdat_file_outliers <- paste0("./Subgroup_MR/sensitivity/Outliers.csv")
outdat_file_distort <- paste0("./Subgroup_MR/sensitivity/MRPRESSO_distort.csv")
write.table(data_mrpresso, file = outdat_file_mrpresso, sep = ",", quote = FALSE, row.names = FALSE)
write.table(data_outliers, file = outdat_file_outliers, sep = ",", quote = FALSE, row.names = FALSE)
write.table(data_distort, file = outdat_file_distort, sep = ",", quote = FALSE, row.names = FALSE)


data <- read.table(paste0("./Subgroup_MR/results_main/MRPRESSO_distort.csv"), 
                   header = T, sep = ",", stringsAsFactors = FALSE)

data_long <- data %>%
  pivot_longer(cols = c(V1, V2), names_to = "variable", values_to = "value")

outdat_file_distort <- paste0("./Subgroup_MR/results_main/MRPRESSO_distort.csv")
write.table(data_long, file = outdat_file_distort, sep = ",", quote = FALSE, row.names = FALSE)


## Repeated MR analysis after removing outliers
data_mo_res <- data.frame()
data_mo_heter <- data.frame()
data_mo_egger <- data.frame()

exposures <- c("LST","MVPA","AbPA","MVPA_2018","ST")
outcomes <- c("ebi-a-GCST006908","ebi-a-GCST006910","ebi-a-GCST006907","ebi-a-GCST006909")

for (exposure in exposures) {
  for (outcome in outcomes) {
    e_o_dat <- read.table(paste0("./Subgroup_MR/data_doutliers/", exposure, "_", outcome, "_V3.csv"), 
                          header = T, sep = ",", stringsAsFactors = FALSE)
    
    print(paste("Step1: Loading:", exposure, " vs ", outcome, "。","No. of SNPs:", nrow(e_o_dat),"."))
    
    e_o_dat_clean <- e_o_dat[e_o_dat$outlier=="no",]
    print(paste("Step2: Dropping outlier SNPs ", nrow(e_o_dat)-nrow(e_o_dat_clean),". Leaving No. of SNPs：",nrow(e_o_dat_clean)))
    
    if (nrow(e_o_dat)-nrow(e_o_dat_clean) > 0) {
      # IVW MR analysis
      res <-generate_odds_ratios(mr(e_o_dat_clean, method_list = c("mr_ivw")))
      # Heterogeneity
      heter <- mr_heterogeneity(e_o_dat_clean)
      # I^2
      heter$I2 <- NA
      heter$I2 <- ifelse(heter$Q < heter$Q_df, 0, paste0(abs(round((heter$Q - heter$Q_df) * 100 / heter$Q, 2)), "%"))
      
      # MR Egger Intercept
      pleio <- mr_pleiotropy_test(e_o_dat_clean)
      
      
      ##Tidy results
      out_mainres <- subset(res, select = c("outcome", "exposure", "method", "nsnp", "b",                                           "se", "pval", "lo_ci", "up_ci", "or", "or_lci95", "or_uci95"))
      out_heter <- subset(heter, select = c("outcome", "exposure", "method", "Q", "Q_df", "Q_pval", "I2"))
      out_pleio <- subset(pleio, select = c("outcome", "exposure", "egger_intercept", "se", "pval"))
      
      exp_name <- switch(exposure,
                         "LST" = "Leisure screen time",
                         "AbPA" = "Accelerometer-based physical activity",
                         "MVPA" = "Leisure time moderate-to-vigorous physical activity",
                         "MVPA_2018" = "Overall moderate-to-vigorous physical activity",
                         "ST" = "Overall sedentary duration",
                         "No-defined")
      out_name <- switch(outcome,
                         "ebi-a-GCST006908" = "Total ischemic stroke",
                         "ebi-a-GCST006910" = "Cardioembolic stroke",
                         "ebi-a-GCST006907" = "Large artery stroke",
                         "ebi-a-GCST006909" = "Small vessel stroke",
                         "No-defined" )
      
      out_heter$outcome <- out_pleio$outcome <- out_name
      out_heter$exposure <- out_pleio$exposure <- exp_name
      
      out_mainres %>%
        mutate(
          heter_p = if_else(method == "Inverse variance weighted", out_heter$Q_pval[out_heter$method == "Inverse variance weighted"], NA_real_),
          eggerint_p = if_else(method == "Inverse variance weighted", out_pleio$pval, NA_real_),
          se = as.numeric(se),
          `Outcomes` = out_name,
          `Exposures` = exp_name,
          `METHOD` = method,
          `SNPs` = nsnp,
          `Beta` = sprintf("%.3f", b),
          `SE` = if_else(is.na(se), "", sprintf("%.5s", se)),
          `OR (95% Cl)` = sprintf("%.2f(%.2f to %.2f)", or, or_lci95, or_uci95),
          `P` = sprintf("%.3f", pval),
          `P_e` = sprintf("%.2e", pval),
          `Q statistic P-value` = if_else(is.na(heter_p), "", sprintf("%.3f", heter_p)),
          `MR-Egger pleiotropy test P-value` = if_else(is.na(eggerint_p), "", sprintf("%.3f", eggerint_p))
        ) -> out_mainres
      
      # Combined results
      data_mo_res <- rbind(data_mo_res, out_mainres)
      data_mo_heter <- rbind(data_mo_heter, out_heter)
      data_mo_egger <- rbind(data_mo_egger, out_pleio)
    } else {
      print(paste("No outliers"))
    }
    
    print(paste("--------------------------------------------------------------------------"))
  }
}


output_file_res <- paste0("./Subgroup_MR/results_moveoutliers/Results_PAtoStroke_mo.csv")
output_file_heter <- paste0("./Subgroup_MR/results_moveoutliers/Heterogeneity_mo.csv")
output_file_pleio <- paste0("./Subgroup_MR/results_moveoutliers/MRegger_intercept_mo.csv")

write.table(data_mo_res, file = output_file_res, sep = ",", quote = FALSE, row.names = FALSE)
write.table(data_mo_heter, file = output_file_heter, sep = ",", quote = FALSE, row.names = FALSE)
write.table(data_mo_egger, file = output_file_pleio, sep = ",", quote = FALSE, row.names = FALSE)
