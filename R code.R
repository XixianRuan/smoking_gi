# The current Mendelian randomization analysis had multiple exposures 
# and outcomes. To avoid duplication codes, the codes only showed the example 
# of analysis on association between smoking initiation and ulcerative colitis.

## Part 1  Mendelian randomization analysis in different datasets ##
library(MRInstruments)
library(readr)
library(TwoSampleMR)

uc_ukb_dat<-read_csv('result_1019/ukb/smoinit/uc.csv')
uc_finn_dat<-read_csv('result_1019/finngen/smoinit/finngen_R7_K11_ULCER.csv')
uc_iibdgc_dat<-read_csv('result_1019/iibdgc/smoinit/ieu-a-32')

# IVW method and heterogeneity test
res_uc_ukb<-mr(uc_ukb_dat, method_list = c('mr_ivw_mre'))
res_uc_finn<-mr(uc_finn_dat, method_list = c('mr_ivw_mre'))
res_uc_iibdgc<-mr(uc_iibdgc_dat, method_list = c('mr_ivw_mre'))

heterogeneity_uc_ukb<-mr_heterogeneity(uc_ukb_dat, method_list = c('mr_ivw_mre'))
heterogeneity_uc_finn<-mr_heterogeneity(uc_finn_dat, method_list = c('mr_ivw_mre'))
heterogeneity_uc_iibdgc<-mr_heterogeneity(uc_iibdgc_dat, method_list = c('mr_ivw_mre'))

# MR-Egger method and pleiotropy test
res_uc_ukb_egger<-mr(uc_ukb_dat, method_list = c('mr_egger_regression'))
res_uc_finn_egger<-mr(uc_finn_dat, method_list = c('mr_egger_regression'))
res_uc_iibdgc_egger<-mr(uc_iibdgc_dat, method_list = c('mr_egger_regression'))

pleiotropy_uc_ukb<-mr_pleiotropy_test(uc_ukb_dat)
pleiotropy_uc_finn<-mr_pleiotropy_test(uc_finn_dat)
pleiotropy_uc_iibdgc<-mr_pleiotropy_test(uc_iibdgc_dat)

# Weighted median method
res_uc_ukb_wm<-mr(uc_ukb_dat, method_list = c('mr_weighted_median'))
res_uc_finn_wm<-mr(uc_finn_dat, method_list = c('mr_weighted_median'))
res_uc_iibdgc_wm<-mr(uc_iibdgc_dat, method_list = c('mr_weighted_median'))

# MR-PRESSO method
presso_uc_ukb<-MRPRESSO::mr_presso(BetaOutcome ='beta.outcome',  BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                             OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = as.data.frame(uc_ukb_dat), NbDistribution = 2000,  
                                             SignifThreshold = 0.05)

presso_uc_finn<-MRPRESSO::mr_presso(BetaOutcome ='beta.outcome',  BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                   OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = as.data.frame(uc_finn_dat), NbDistribution = 2000,  
                                   SignifThreshold = 0.05)

presso_uc_iibdgc<-MRPRESSO::mr_presso(BetaOutcome ='beta.outcome',  BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                   OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = as.data.frame(uc_iibdgc_dat), NbDistribution = 2000,  
                                   SignifThreshold = 0.05)


## Part 2 Meta-analysis by combining estimates from different datasets
library(metafor)
uc_ivw_meta<-rma(yi=b,
                 sei = se,
                 data=rbind(res_uc_ukb, res_uc_finn, res_uc_iibdgc),
                 method = 'FE')
summary(uc_ivw_meta)

uc_egger_meta<-rma(yi=b,
                 sei = se,
                 data=rbind(res_uc_ukb_egger, res_uc_finn_egger, res_uc_iibdgc_egger),
                 method = 'FE')
summary(uc_egger_meta)

uc_wm_meta<-rma(yi=b,
                   sei = se,
                   data=rbind(res_uc_ukb_wm, res_uc_finn_wm, res_uc_iibdgc_wm),
                   method = 'FE')
summary(uc_wm_meta)

