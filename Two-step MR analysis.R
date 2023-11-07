########################### Two step MR #########################
########################### Date: 20231021 ####################
############################ Author: Yidie #####################

#Body mass index ieu-b-40
#Waist to hip ratio 
#Body fat percentage ukb-b-8909
#Waist circumference ukb-b-9405
#Hip circumference ukb-b-15590


########################################First step######################################
###########################Part1 Create datasets for analysis########################

library(TwoSampleMR)
library(plyr)
library(MRPRESSO)
library(mr.raps)
library(tidyverse)

##Part one--Prepare data
exposures <- c("LST","MVPA","AbPA")

##obesity and lipids
outcomes <- c("ieu-b-40","ukb-b-8909","ukb-b-9405","ukb-b-15590","ieu-b-107","ieu-b-108")


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
    
    ##step 5: Drop ambiguous SNPs
    cleaned_data <- e_o_dat_clean[e_o_dat_clean$ambiguous != "TRUE", ]

    #Output harmonised data
    outdat_file <- paste0("./Twostep_MR/First_step/data_V1/data_origin/", exposure1, "_", outcome1, "_V1.csv")
    write.table(cleaned_data, file = outdat_file, sep = ",", quote = FALSE, row.names = FALSE)
    ##
    ##
    print("END: e_o_dat_clean prepared for two-sample MR analysis")
    
    print("--------------------------------------------------------------------------")
    
  }
  
}




exposures <- c("LST","MVPA","AbPA")
outcomes <- c("WHR","LDL_C","HDL_C","TG")
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
    outcome_dat <- read.table(paste0("./match_V3/", exposure1, "_", outcome1, ".matchV3.tsv"), 
                              header = T, sep = "\t", stringsAsFactors = FALSE)
    
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
    
    ##step 5: Drop ambiguous SNPs
    cleaned_data <- e_o_dat_clean[e_o_dat_clean$ambiguous != "TRUE", ]
    
    #输出Harmonised数据
    outdat_file <- paste0("./Twostep_MR/First_step/data_V1/data_origin/", exposure1, "_", outcome1, "_V1.csv")
    write.table(cleaned_data, file = outdat_file, sep = ",", quote = FALSE, row.names = FALSE)
    ##
    print("END: cleaned_data prepared for two-sample MR analysis")
    
    print("--------------------------------------------------------------------------")
    
  }
}



###########################Part 2. MR from exposures to mediatos ########################

data_res <- data.frame()
data_heter <- data.frame()
data_egger <- data.frame()
#exposures and outcomes
exposures <- c("MVPA","AbPA","LST")
outcomes <- c("WHR","ieu-b-40","ukb-b-8909","ukb-b-9405","ukb-b-15590","ieu-b-107","ieu-b-108","LDL_C","HDL_C","TG")

#### Primary analysis
for (exposure in exposures) {
  for (outcome in outcomes) {
    # read formated dataset
    mydat <- read.delim(file = paste0("./Twostep_MR/First_step/data_V1/data_origin/", exposure, "_", outcome, "_V1.csv"), header = TRUE, sep = ",")
    
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
    print(paste0("step3 mr-raps: DONE"))
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
    out_name <- switch(outcome1,
                       "WHR" = "Waist to hip ratio",
                       "ieu-b-40" = "Body mass index",
                       "ukb-b-8909" = "Body fat percentage",
                       "ukb-b-9405" = "Waist circumference",
                       "ukb-b-15590" = "Hip circumference",
                       "LDL_C" = "Low density lipoprotein cholesterol",
                       "HDL_C" = "High density lipoprotein cholesterol",
                       "ieu-b-107" = "Apolipoprotein A-1",
                       "ieu-b-108" = "Apolipoprotein B",
                       "TG" = "Triglyceride levels",
                       "No-defined"  
    )
    
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


output_file_res <- paste0("./Twostep_MR/First_step/results_V1/Results_PAtoMediatorsres.csv")
output_file_heter <- paste0("./Twostep_MR/First_step/results_V1/Heterogeneity.csv")
output_file_pleio <- paste0("./Twostep_MR/First_step/results_V1/MRegger_intercept.csv")

write.table(data_res, file = output_file_res, sep = ",", quote = FALSE, row.names = FALSE)
write.table(data_heter, file = output_file_heter, sep = ",", quote = FALSE, row.names = FALSE)
write.table(data_egger, file = output_file_pleio, sep = ",", quote = FALSE, row.names = FALSE)






#####################################################################
##################Part 3. Detecting outliers#######################
data_mrpresso <- data.frame()
data_outliers <- data.frame()
data_distort <- data.frame()

#exposures and outcomes
exposures <- c("LST","MVPA","AbPA")
outcomes <- c("WHR","ieu-b-40","ukb-b-8909","ukb-b-9405","ukb-b-15590","ieu-b-107","ieu-b-108","LDL_C","HDL_C","TG")

#MR PRESSO test
for (exposure1 in exposures) {
  for (outcome1 in outcomes) {
    #step 1: prepare the exposure data
    e_o_dat_clean <- read.table(paste0("./Twostep_MR/First_step/data_V1/data_origin/", exposure1, "_", outcome1, "_V1.csv"), 
                                header = T, sep = ",", stringsAsFactors = FALSE)
    print(paste0("Loading datasets：", exposure1,"_", outcome1,". Time：", Sys.time ()))
    
    exp_name <- switch(exposure1,
                       "LST" = "Leisure screen time",
                       "AbPA" = "Accelerometer-based physical activity",
                       "MVPA" = "Leisure time moderate-to-vigorous physical activity",
                       "No-defined"  
    )
    out_name <- switch(outcome1,
                       "WHR" = "Waist to hip ratio",
                       "ieu-b-40" = "Body mass index",
                       "ukb-b-8909" = "Body fat percentage",
                       "ukb-b-9405" = "Waist circumference",
                       "ukb-b-15590" = "Hip circumference",
                       "LDL_C" = "Low density lipoprotein cholesterol",
                       "HDL_C" = "High density lipoprotein cholesterol",
                       "ieu-b-107" = "Apolipoprotein A-1",
                       "ieu-b-108" = "Apolipoprotein B",
                       "TG" = "Triglyceride levels",
                       "No-defined"  
    )
    
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
    outdat_file <- paste0("./Twostep_MR/First_step/data_V1/data_mrpresso/", exposure1, "_", outcome1, "_mrpresso.csv")
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
outdat_file_mrpresso <- paste0("./Twostep_MR/First_step/results_V1/MRPRESSO_lipids.csv")
outdat_file_outliers <- paste0("./Twostep_MR/First_step/results_V1/Outliers_lipids.csv")
outdat_file_distort <- paste0("./Twostep_MR/First_step/results_V1/MRPRESSO_distort_lipids.csv")
write.table(data_mrpresso, file = outdat_file_mrpresso, sep = ",", quote = FALSE, row.names = FALSE)
write.table(data_outliers, file = outdat_file_outliers, sep = ",", quote = FALSE, row.names = FALSE)
write.table(data_distort, file = outdat_file_distort, sep = ",", quote = FALSE, row.names = FALSE)


data <- read.table(paste0("./Twostep_MR/First_step/results_V1/MRPRESSO_distort_lipids.csv"), 
                   header = T, sep = ",", stringsAsFactors = FALSE)
data$V1 <- as.character(data$V1)
data_long <- data %>%
  pivot_longer(cols = c(V1, V2), names_to = "variable", values_to = "value")

outdat_file_distort <- paste0("./Twostep_MR/First_step/results_V1/MRPRESSO_distort_lipids.csv")
write.table(data_long, file = outdat_file_distort, sep = ",", quote = FALSE, row.names = FALSE)

data <- read.table(paste0("./Twostep_MR/First_step/results_V1/MRPRESSO_distort_obesity.csv"), 
                   header = T, sep = ",", stringsAsFactors = FALSE)

data_long <- data %>%
  pivot_longer(cols = c(V1, V2), names_to = "variable", values_to = "value")

outdat_file_distort <- paste0("./Twostep_MR/First_step/results_V1/MRPRESSO_distort_obesity.csv")
write.table(data_long, file = outdat_file_distort, sep = ",", quote = FALSE, row.names = FALSE)




## Repeated MR analysis after removing outliers
data_mo_res <- data.frame()
data_mo_heter <- data.frame()
data_mo_egger <- data.frame()

exposures <- c("LST","MVPA","AbPA")
outcomes <- c("WHR","ieu-b-40","ukb-b-8909","ukb-b-9405","ukb-b-15590","ieu-b-107","ieu-b-108","LDL_C","HDL_C","TG")

for (exposure in exposures) {
  for (outcome in outcomes) {
    
    try(e_o_dat <- read.table(paste0("./Twostep_MR/First_step/data_V1/data_mrpresso/", exposure, "_", outcome, "_mrpresso.csv"), 
                              header = T, sep = ",", stringsAsFactors = FALSE), silent = FALSE,
        outFile = getOption("try.outFile", default = stderr()))
   
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
                         "No-defined")
      out_name <- switch(outcome1,
                         "WHR" = "Waist to hip ratio",
                         "ieu-b-40" = "Body mass index",
                         "ukb-b-8909" = "Body fat percentage",
                         "ukb-b-9405" = "Waist circumference",
                         "ukb-b-15590" = "Hip circumference",
                         "LDL_C" = "Low density lipoprotein cholesterol",
                         "HDL_C" = "High density lipoprotein cholesterol",
                         "ieu-b-107" = "Apolipoprotein A-1",
                         "ieu-b-108" = "Apolipoprotein B",
                         "TG" = "Triglyceride levels",
                         "No-defined"  
      )
      
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

output_file_res <- paste0("./Twostep_MR/First_step/results_V1/Results_PAtoMediators_mo.csv")
output_file_heter <- paste0("./Twostep_MR/First_step/results_V1/Heterogeneity_mo.csv")
output_file_pleio <- paste0("./Twostep_MR/First_step/results_V1/MRegger_intercept_mo.csv")

write.table(data_mo_res, file = output_file_res, sep = ",", quote = FALSE, row.names = FALSE)
write.table(data_mo_heter, file = output_file_heter, sep = ",", quote = FALSE, row.names = FALSE)
write.table(data_mo_egger, file = output_file_pleio, sep = ",", quote = FALSE, row.names = FALSE)




###Output the SNPs lists that used in STEP ONE
d_list <- data.frame()
exposures <- c("LST","MVPA","AbPA")
outcomes <- c("WHR","ieu-b-40","ukb-b-8909","ukb-b-9405","ukb-b-15590","ieu-b-107","ieu-b-108","LDL_C","HDL_C","TG")

for (exposure in exposures) {
  for (outcome in outcomes) {
    path <- paste0("./Twostep_MR/First_step/data_doutliers/", exposure, "_", outcome, "_V2.csv")
    d <- read.csv(file=path, header = TRUE)
    d <- as.data.frame(d[, c("SNP")])
    d_list <- rbind(d_list, d)
    
  }
}
colnames(d_list)[1] <- 'SNP'
output_file <- paste0("./Twostep_MR/",outcome,".txt")

write.table(d_list, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)

rm("d_list", "path", "d", "exposure","exposures","outcome","outcomes","output_file")









########################################Second step######################################
###########################Part1 Create datasets for analysis########################


##Part one--Prepare data
#exposures and outcomes
exposures <- c("BMI","BFP","WHR","WC","HC","LDL_C","HDL_C","TG","ApA1","ApB")
outcomes <- c("AIS","LAS","CAS","SVS")
for (exposure1 in exposures) {
  for (outcome1 in outcomes) {
    
    print(paste0("BEGIN: prepare data ", exposure1,"_", outcome1))
    #prepare data
    #step 1: prepare the exposure data
    exposure <- read.table(paste0("./clumped_plink/",exposure1,".clumped.tsv"), 
                           header = T, sep = "\t", stringsAsFactors = FALSE)
    
    #orienta to the positive direction
    exposure$EA_c <- exposure$EA; exposure$NEA_c <- exposure$NEA
    exposure[ which(exposure$BETA < 0), "EA"] <-   exposure[ which(exposure$BETA < 0), "NEA_c"] 
    exposure[ which(exposure$BETA < 0), "NEA"] <-   exposure[ which(exposure$BETA < 0), "EA_c"]  
    exposure[ which(exposure$BETA < 0), "EAF"] <-  1 - exposure[ which(exposure$BETA < 0), "EAF"]  
    exposure$BETA <- abs(exposure$BETA)  
    exposure <- exposure[, c(1:11)] # SNP CHR POS EA NEA EAF BETA SE P. need to check the EAF colunm before
    ##
    exposure_dat <- format_data(exposure, type = "exposure", header = TRUE, phenotype_col = "PHENOTYPE",
                                snp_col = "SNP", beta_col = "BETA", se_col = "SE", eaf_col = "EAF", effect_allele_col = "EA",
                                other_allele_col = "NEA", pval_col = "P", samplesize_col = "N")
    exposure_dat_tmp <- exposure_dat
    print(paste("step 1. ", nrow(exposure_dat_tmp), "SNPs", ". Data dimension: ", dim(exposure_dat_tmp)[1], dim(exposure_dat_tmp)[2]))
    ##
    
    #step 2: prepare the outcome data by extract the summary statistics using above IVs information
    outcome_dat <- read.table(paste0("./match_V3/", exposure1, "_", outcome1, ".matchV3.tsv"), 
                              header = T, sep = "\t", stringsAsFactors = FALSE)
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
    
    ##Dropping ambigurous SNPs
    cleaned_data <- e_o_dat_clean[e_o_dat_clean$ambiguous != "TRUE", ]
    print(paste("step 4. Drop ", exposure1,"_", outcome1, "ambigurous SNPs：", e_o_dat_clean$SNP[e_o_dat_clean$ambiguous == "TRUE"], ". Leaving：", nrow(cleaned_data), " SNPs."))
    
    ##step5 Removeing SNPs used as IVs in the First step MR
    drop_list <- read.table(paste0("./Twostep_MR/", exposure1, ".txt"), 
                            header = T, sep = "\t", stringsAsFactors = FALSE)
    
    
    cleaned_data <- cleaned_data[!(cleaned_data$SNP %in% drop_list$SNP), ]
    print(paste("step 4. Removeing SNPs used as IVs in the First step MR. Leaving: ", nrow(cleaned_data), "."))
    ##
    if (exposure1 == "ApA1") {
      cleaned_data$samplesize.exposure <- 393193
    } else if (exposure1 == "ApB") {
      cleaned_data$samplesize.exposure <- 439214
    }
    cleaned_data2 <- steiger_filtering(cleaned_data)
    #Output Harmonised data
    outdat_file <- paste0("./Twostep_MR/Second_step/data/", exposure1, "_", outcome1, "_V1.csv")
    write.table(cleaned_data2, file = outdat_file, sep = ",", quote = FALSE, row.names = FALSE)
    ##
    print("END: cleaned_data2 prepared for two-sample MR analysis")
    
    print("-------------------------------------------------------------------------")
    
  }
  
}



###########################Part 2. MR-From mediators to outcomes########################

data_res <- data.frame()
data_heter <- data.frame()
data_egger <- data.frame()

#exposures and outcomes
exposures <- c("BMI","BFP","WHR","WC","HC","LDL_C","HDL_C","TG","ApA1","ApB")
outcomes <- c("AIS","LAS","CAS","SVS")

#### Primary analysis
for (exposure in exposures) {
  for (outcome in outcomes) {
    # 读取和格式化数据
    mydat <- read.delim(file = paste0("./Twostep_MR/Second_step/data/", exposure, "_", outcome, "_V1.csv"), header = TRUE, sep = ",")
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
    
    exp_name <- switch(exposure1,
                       "WHR" = "Waist to hip ratio",
                       "ieu-b-40" = "Body mass index",
                       "ukb-b-8909" = "Body fat percentage",
                       "ukb-b-9405" = "Waist circumference",
                       "ukb-b-15590" = "Hip circumference",
                       "LDL_C" = "Low density lipoprotein cholesterol",
                       "HDL_C" = "High density lipoprotein cholesterol",
                       "ieu-b-107" = "Apolipoprotein A-1",
                       "ieu-b-108" = "Apolipoprotein B",
                       "TG" = "Triglyceride levels",
                       "No-defined"  
    )
    out_name <- switch(outcome,
                       "AIS" = "Total ischemic stroke",
                       "CAS" = "Cardioembolic stroke",
                       "LAS" = "Large artery stroke",
                       "SVS" = "Small vessel stroke",
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

output_file_res <- paste0("./Twostep_MR/Second_step/results/Results_Mediators_to_Stroke.csv")
output_file_heter <- paste0("./Twostep_MR/Second_step/results/Heterogeneity.csv")
output_file_pleio <- paste0("./Twostep_MR/Second_step/results/MRegger_intercept.csv")

write.table(data_res, file = output_file_res, sep = ",", quote = FALSE, row.names = FALSE)
write.table(data_heter, file = output_file_heter, sep = ",", quote = FALSE, row.names = FALSE)
write.table(data_egger, file = output_file_pleio, sep = ",", quote = FALSE, row.names = FALSE)



#####################################################################
##################Part 3. 检验是否存在outliers#######################
data_mrpresso <- data.frame()
data_outliers <- data.frame()
data_distort <- data.frame()
#exposures and outcomes
exposures <- c("BMI","BFP","WHR","WC","HC","LDL_C","HDL_C","TG","ApA1","ApB")
outcomes <- c("AIS","LAS","CAS","SVS")

#MR PRESSO test
for (exposure1 in exposures) {
  for (outcome1 in outcomes) {
    #step 1: prepare the exposure data
    e_o_dat_clean <- read.table(paste0("./Twostep_MR/Second_step/data/", exposure1, "_", outcome1, "_V1.csv"), 
                                header = T, sep = ",", stringsAsFactors = FALSE)
    print(paste0("Loading datasets：", exposure1,"_", outcome1,". Time：", Sys.time ()))
    exp_name <- switch(exposure1,
                       "WHR" = "Waist to hip ratio",
                       "BMI" = "Body mass index",
                       "BFP" = "Body fat percentage",
                       "WC" = "Waist circumference",
                       "HC" = "Hip circumference",
                       "LDL_C" = "Low density lipoprotein cholesterol",
                       "HDL_C" = "High density lipoprotein cholesterol",
                       "ApA1" = "Apolipoprotein A-1",
                       "ApB" = "Apolipoprotein B",
                       "TG" = "Triglyceride levels",
                       "No-defined"  
    )
    out_name <- switch(outcome1,
                       "AIS" = "Total ischemic stroke",
                       "CAS" = "Cardioembolic stroke",
                       "LAS" = "Large artery stroke",
                       "SVS" = "Small vessel stroke",
                       "No-defined" 
    )
    
    
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
    outdat_file <- paste0("./Twostep_MR/Second_step/data_doutliers/", exposure1, "_", outcome1, "_V2.csv")
    write.table(e_o_dat_clean, file = outdat_file, sep = ",", quote = FALSE, row.names = FALSE)
    ##
    print(paste("END: MRPRESSO detected ", exposure1, "-", outcome1," outliers. ", e_o_dat_clean$SNP[e_o_dat_clean$outlier == "yes"]))
    print(paste("---------------------------------------------------------------------------"))
    print(paste("DONE：", Sys.time ()))
    ##
    gc()
  }
  
}
colnames(data_mrpresso) <- c("exposue","outcome","MR Analysis", "Causal Estimate", "Sd",  "T-stat","P-value","RSSobs")

data_mrpresso <- cbind(data_mrpresso[1:7], unlist(data_mrpresso$RSSobs))
outdat_file_mrpresso <- paste0("./Twostep_MR/Second_step/results/MRPRESSO.csv")
outdat_file_outliers <- paste0("./Twostep_MR/Second_step/results/Outliers.csv")
outdat_file_distort <- paste0("./Twostep_MR/Second_step/results/MRPRESSO_distort.csv")
write.table(data_mrpresso, file = outdat_file_mrpresso, sep = ",", quote = FALSE, row.names = FALSE)
write.table(data_outliers, file = outdat_file_outliers, sep = ",", quote = FALSE, row.names = FALSE)
write.table(data_distort, file = outdat_file_distort, sep = ",", quote = FALSE, row.names = FALSE)


data <- read.table(paste0("./Twostep_MR/Second_step/results/MRPRESSO_distort.csv"), 
                   header = T, sep = ",", stringsAsFactors = FALSE)

data_long <- data %>%
  pivot_longer(cols = c(V1, V2), names_to = "variable", values_to = "value")

outdat_file_distort <- paste0("./Twostep_MR/Second_step/results/MRPRESSO_distort.csv")
write.table(data_long, file = outdat_file_distort, sep = ",", quote = FALSE, row.names = FALSE)





## Repeated MR analysis after removing outliers
data_mo_res <- data.frame()
data_mo_heter <- data.frame()
data_mo_egger <- data.frame()

exposures <- c("BMI","BFP","WHR","WC","HC","LDL_C","HDL_C","TG","ApA1","ApB")
outcomes <- c("AIS","LAS","CAS","SVS")
for (exposure in exposures) {
  for (outcome in outcomes) {
    e_o_dat <- read.table(paste0("./Twostep_MR/Second_step/data_doutliers/", exposure, "_", outcome, "_V2.csv"), 
                          header = T, sep = ",", stringsAsFactors = FALSE)
    
    print(paste("Step1: Loading:", exposure, " vs ", outcome, "。","No. of SNPs:", nrow(e_o_dat),"."))
    
    e_o_dat_clean <- e_o_dat[e_o_dat$outlier=="no",]
    print(paste("Step2: Dropping outlier SNPs ", nrow(e_o_dat)-nrow(e_o_dat_clean),". Leaving No. of SNPs：",nrow(e_o_dat_clean)))
    
    
    if (nrow(e_o_dat)-nrow(e_o_dat_clean) > 0) {
      ## removing ambigurous SNPs
      cleaned_data <- e_o_dat_clean[e_o_dat_clean$ambiguous != "TRUE", ]
      ## removing SNPs used in the First step
      drop_list <- read.table(paste0("./Twostep_MR/", exposure, ".txt"), 
                              header = T, sep = "\t", stringsAsFactors = FALSE)
      
      
      cleaned_data <- cleaned_data[!(cleaned_data$SNP %in% drop_list$SNP), ]
      
      cleaned_data <- cleaned_data %>%filter(pval.outcome > 5e-06)
      
      e_o_dat_clean <- cleaned_data
      print(paste("e_o_dat_clean has：", nrow(e_o_dat_clean), " SNPs."))
      
      
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
                         "WHR" = "Waist to hip ratio",
                         "BMI" = "Body mass index",
                         "BFP" = "Body fat percentage",
                         "WC" = "Waist circumference",
                         "HC" = "Hip circumference",
                         "LDL_C" = "Low density lipoprotein cholesterol",
                         "HDL_C" = "High density lipoprotein cholesterol",
                         "ApA1" = "Apolipoprotein A-1",
                         "ApB" = "Apolipoprotein B",
                         "TG" = "Triglyceride levels",
                         "No-defined")
      out_name <- switch(outcome,
                         "AIS" = "Total ischemic stroke",
                         "CAS" = "Cardioembolic stroke",
                         "LAS" = "Large artery stroke",
                         "SVS" = "Small vessel stroke",
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
    
    print(paste("-------------------------------END：分割线--------------------------------------"))
  }
}

output_file_res <- paste0("./Twostep_MR/Second_step/results_moveoutliers/Results_Mediators_to_Stroke_mo.csv")
output_file_heter <- paste0("./Twostep_MR/Second_step/results_moveoutliers/Heterogeneity_mo.csv")
output_file_pleio <- paste0("./Twostep_MR/Second_step/results_moveoutliers/MRegger_intercept_mo.csv")

write.table(data_mo_res, file = output_file_res, sep = ",", quote = FALSE, row.names = FALSE)
write.table(data_mo_heter, file = output_file_heter, sep = ",", quote = FALSE, row.names = FALSE)
write.table(data_mo_egger, file = output_file_pleio, sep = ",", quote = FALSE, row.names = FALSE)








############################# Mediation analysis #################################

#.rs.restartR()

##Obesity related
data <- data.frame()
exposures <- c("LST","MVPA","AbPA")
mediators <- c("BMI","BFP","WHR","WC","HC")
outcomes <- c("AIS","CAS","LAS","SVS")

for (exposure in exposures) {
  for (mediator in mediators) {
    for (outcome in outcomes) {
      
      print(paste0("BEGIN: -------prepare data ", exposure,"_", mediator,"_",outcome,"---------"))
      
      if (mediator == "BMI") {
        step1_res <- read.csv(paste0("./Twostep_MR/First_step/results_V1/mrpresso/", exposure, "_ieu-b-40_mrpresso.csv"), 
                              header = T, sep = ",", stringsAsFactors = FALSE)
      } else if (mediator == "BFP"){
        step1_res <- read.csv(paste0("./Twostep_MR/First_step/results_V1/mrpresso/", exposure, "_ukb-b-8909_mrpresso.csv"), 
                              header = T, sep = ",", stringsAsFactors = FALSE)
      } else if (mediator == "WHR"){
        step1_res <- read.csv(paste0("./Twostep_MR/First_step/results_V1/mrpresso/", exposure, "_WHR_mrpresso.csv"), 
                              header = T, sep = ",", stringsAsFactors = FALSE)
      } else if (mediator == "WC"){
        step1_res <- read.csv(paste0("./Twostep_MR/First_step/results_V1/mrpresso/", exposure, "_ukb-b-9405_mrpresso.csv"), 
                              header = T, sep = ",", stringsAsFactors = FALSE)
      } else if (mediator == "HC"){
        step1_res <- read.csv(paste0("./Twostep_MR/First_step/results_V1/mrpresso/", exposure, "_ukb-b-15590_mrpresso.csv"), 
                              header = T, sep = ",", stringsAsFactors = FALSE)
      }
      
      step2_res <- read.csv(paste0("./Twostep_MR/Second_step/results_moveoutliers/", mediator, "_", outcome, "_mrpresso.csv"), 
                            header = T, sep = ",", stringsAsFactors = FALSE)
      
      res_2smr <- read.csv(paste0("./Forward_MR/results_moveoutliers/", exposure, "_", outcome, "_mrpresso.csv"), 
                           header = T, sep = ",", stringsAsFactors = FALSE)
      
      step1res <- step1_res[step1_res$method=="Inverse variance weighted",]
      step2res <- step2_res[step2_res$method=="Inverse variance weighted",] 
      res2smr <- res_2smr[res_2smr$method=="Inverse variance weighted",] 
      
      d_res <- subset(step1res, select = c("exposure"))
      d_res$mediator <- step1res$outcome
      d_res$outcome <- step2res$outcome
      d_res$nsnp1 <- step1res$nsnp
      d_res$nsnp2 <- step2res$nsnp
      d_res$b1 <- step1res$b
      d_res$b2 <- step2res$b
      d_res$se1 <- step1res$se
      d_res$se2 <- step2res$se
      
      #####data prepared#########
      print(paste0("step1: data prepared ", exposure,"_", mediator,"_",outcome,"。"))
      
      
      ####Calculating indirect effect and its se #####
      d_res$b3 <- d_res$b1*d_res$b2
      d_res$se3 <- sqrt(d_res$b1*d_res$b1*d_res$se2*d_res$se2+d_res$b2*d_res$b2*d_res$se1*d_res$se1)
      confidence_level <- 0.95
      margin_of_error <- qnorm((1 + confidence_level) / 2) * d_res$se3
      d_res$b3_LL <- d_res$b3 - margin_of_error
      d_res$b3_UL <- d_res$b3 + margin_of_error
      
      d_res$Z3 <- d_res$b3/d_res$se3
      d_res$p3 <-  2*pnorm(q=abs(d_res$Z3), lower.tail=FALSE)
      
      print(paste0("step2: indirect effect is ",d_res$b3,", p value is ",d_res$p3,"."))
      
      #### Propartion mediated #####
      d_res$mp = 100*d_res$b3/res2smr$b
      d_res$mp_LL <- 100*d_res$b3_LL/res2smr$b
      d_res$mp_UL <- 100*d_res$b3_UL/res2smr$b
      
      ####Tidy results##
      
      d_res$total <- sprintf("%.3f(%.3f to %.3f)",
                             res2smr$b,
                             res2smr$b-qnorm((1 + 0.95) / 2) * res2smr$se,
                             res2smr$b+qnorm((1 + 0.95) / 2) * res2smr$se)
      
      d_res$step1 <- sprintf("%.3f(%.3f to %.3f)",
                             d_res$b1,
                             d_res$b1-qnorm((1 + 0.95) / 2) * d_res$se1,
                             d_res$b1+qnorm((1 + 0.95) / 2) * d_res$se1)
      
      d_res$step2 <- sprintf("%.3f(%.3f to %.3f)",
                             d_res$b2,
                             d_res$b2-qnorm((1 + 0.95) / 2) * d_res$se2,
                             d_res$b2+qnorm((1 + 0.95) / 2) * d_res$se2)
      
      d_res$indirect <- sprintf("%.3f(%.3f to %.3f)",
                                d_res$b3,
                                d_res$b3_LL,
                                d_res$b3_UL)
      
      if (exposure == "MVPA"){
        d_res$Proportion_mediated <- sprintf("%.2f%%(%.2f%% to %.2f%%)", 
                                             d_res$mp,
                                             d_res$mp_UL,
                                             d_res$mp_LL)
      } else if (exposure == "LST"){
        d_res$Proportion_mediated <- sprintf("%.2f%%(%.2f%% to %.2f%%)", 
                                             d_res$mp,
                                             d_res$mp_LL,
                                             d_res$mp_UL)  
      }
      print(paste0("step3: proportion mediated is ",d_res$Proportion_mediated,"。"))
      
      
      data <- rbind(data, d_res)
      
      print(paste0("step4: results combined。"))
      print(paste0("END: ------------------------------------"))
      
    }
  }
}

output_file <- paste0("./Twostep_MR/mediation_results_obesity.csv")

write.table(data, file = output_file, sep = ",", quote = FALSE, row.names = FALSE)
#d_res$seMP <- (d_res$se3^2+d_res$Proportion_mediated^2*res2smr$se^2)/res2smr$b^2


###Lipids
data <- data.frame()
exposures <- c("LST","MVPA","AbPA")
mediators <- c("ApA1","ApB","LDL_C","HDL_C","TG")
outcomes <- c("AIS","CAS","LAS","SVS")

for (exposure in exposures) {
  for (mediator in mediators) {
    for (outcome in outcomes) {
      
      print(paste0("BEGIN: -------prepare data ", exposure,"_", mediator,"_",outcome,"---------"))
      
      if (mediator == "ApA1") {
        step1_res <- read.csv(paste0("./Twostep_MR/First_step/results_V1/mrpresso/", exposure, "_ieu-b-107_mrpresso.csv"), 
                              header = T, sep = ",", stringsAsFactors = FALSE)
      } else if (mediator == "ApB"){
        step1_res <- read.csv(paste0("./Twostep_MR/First_step/results_V1/mrpresso/", exposure, "_ieu-b-108_mrpresso.csv"), 
                              header = T, sep = ",", stringsAsFactors = FALSE)
      } else {
        step1_res <- read.csv(paste0("./Twostep_MR/First_step/results_V1/mrpresso/", exposure, "_",mediator,"_mrpresso.csv"), 
                              header = T, sep = ",", stringsAsFactors = FALSE)
      } 
      
      step2_res <- read.csv(paste0("./Twostep_MR/Second_step/results_moveoutliers/", mediator, "_", outcome, "_mrpresso.csv"), 
                            header = T, sep = ",", stringsAsFactors = FALSE)
      
      res_2smr <- read.csv(paste0("./Forward_MR/results_moveoutliers/", exposure, "_", outcome, "_mrpresso.csv"), 
                           header = T, sep = ",", stringsAsFactors = FALSE)
      
      step1res <- step1_res[step1_res$method=="Inverse variance weighted",]
      step2res <- step2_res[step2_res$method=="Inverse variance weighted",] 
      res2smr <- res_2smr[res_2smr$method=="Inverse variance weighted",] 
      
      d_res <- subset(step1res, select = c("exposure"))
      d_res$mediator <- step1res$outcome
      d_res$outcome <- step2res$outcome
      d_res$nsnp1 <- step1res$nsnp
      d_res$nsnp2 <- step2res$nsnp
      d_res$b1 <- step1res$b
      d_res$b2 <- step2res$b
      d_res$se1 <- step1res$se
      d_res$se2 <- step2res$se
      
      #####Data prepared#########
      print(paste0("step1: data prepared ", exposure,"_", mediator,"_",outcome,"。"))
      
      
      ####Calculating indirect effect and its se #####
      d_res$b3 <- d_res$b1*d_res$b2
      d_res$se3 <- sqrt(d_res$b1*d_res$b1*d_res$se2*d_res$se2+d_res$b2*d_res$b2*d_res$se1*d_res$se1)
      confidence_level <- 0.95
      margin_of_error <- qnorm((1 + confidence_level) / 2) * d_res$se3
      d_res$b3_LL <- d_res$b3 - margin_of_error
      d_res$b3_UL <- d_res$b3 + margin_of_error
      
      d_res$Z3 <- d_res$b3/d_res$se3
      d_res$p3 <-  2*pnorm(q=abs(d_res$Z3), lower.tail=FALSE)
      
      print(paste0("step2: indirect effect is",d_res$b3,". p value is ",d_res$p3,"."))
      
      ####Proportion mediated #####
      d_res$mp = 100*d_res$b3/res2smr$b
      d_res$mp_LL <- 100*d_res$b3_LL/res2smr$b
      d_res$mp_UL <- 100*d_res$b3_UL/res2smr$b
      
      ####Tidy results###
      
      d_res$total <- sprintf("%.3f(%.3f to %.3f)",
                             res2smr$b,
                             res2smr$b-qnorm((1 + 0.95) / 2) * res2smr$se,
                             res2smr$b+qnorm((1 + 0.95) / 2) * res2smr$se)
      
      d_res$step1 <- sprintf("%.3f(%.3f to %.3f)",
                             d_res$b1,
                             d_res$b1-qnorm((1 + 0.95) / 2) * d_res$se1,
                             d_res$b1+qnorm((1 + 0.95) / 2) * d_res$se1)
      
      d_res$step2 <- sprintf("%.3f(%.3f to %.3f)",
                             d_res$b2,
                             d_res$b2-qnorm((1 + 0.95) / 2) * d_res$se2,
                             d_res$b2+qnorm((1 + 0.95) / 2) * d_res$se2)
      
      d_res$indirect <- sprintf("%.3f(%.3f to %.3f)",
                                d_res$b3,
                                d_res$b3_LL,
                                d_res$b3_UL)
      
      if (exposure == "MVPA" ){
        d_res$Proportion_mediated <- sprintf("%.2f%%(%.2f%% to %.2f%%)", 
                                             d_res$mp,
                                             d_res$mp_UL,
                                             d_res$mp_LL)
      } else if (exposure == "LST"){
        d_res$Proportion_mediated <- sprintf("%.2f%%(%.2f%% to %.2f%%)", 
                                             d_res$mp,
                                             d_res$mp_LL,
                                             d_res$mp_UL)  
      } else if (exposure == "AbPA"){
        d_res$Proportion_mediated <- sprintf("%.2f%%(%.2f%% to %.2f%%)", 
                                             d_res$mp,
                                             d_res$mp_UL,
                                             d_res$mp_LL)  
      }
      print(paste0("step3: proportion mediated is ",d_res$Proportion_mediated,"。"))
      
      
      data <- rbind(data, d_res)
      
      print(paste0("step4: results combined。"))
      print(paste0("END: ----------------------------------"))
      
    }
  }
}

output_file <- paste0("./Twostep_MR/mediation_results_lipids.csv")

write.table(data, file = output_file, sep = ",", quote = FALSE, row.names = FALSE)





