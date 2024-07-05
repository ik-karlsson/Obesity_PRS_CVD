################################################################################################
# Name: PRS_BMI_posthocQC_20230123
# Purpose: Compare the predictive value of the PC adjusted PRS to the unadjusted + PCs as cov
#
# Study: SfoEpi: GPdisc_CVD
#
# Created by Ida Karlsson (IK)
# Department of Medical Epidemiology and Biostatistics (MEB), Karolinska Institutet, Stockholm
#
# Created: 	20230123 by Ida Karlsson (IK)
# Updated: 	
#
################################################################################################

.libPaths("Z:/Programs/R/Packages")
# Package to read sas files
library(haven)
# Package for getting coefficients
library(lmtest)
# Package to adjust analyses for clustered data
library(multiwayvcov)

##### Import PCs
pc <- read.table("P:/Dementia_IK/PC_TG_merge/PCs_GHOST/GHOSTY/PCs_GHOSTY_20211005.txt", header=T, sep="\t")
# First 5
pc <- pc[c(1:6,12)]

########## GOSH/PsychChip
##### BMI PRS
BMIprs_GOSH <- read.table("P:/IGEMS/IGEMS_Research/PRS_work/Data/BMI_Yengo2018_SBayesR_STR_GOSH.txt", header=T, sep="\t")
##### Standardize and add study variable
BMIprs_GOSH$zPRS_BMI <- scale(BMIprs_GOSH$BMI_Yengo2018)
BMIprs_GOSH$study <- 'GOSH'

########## TwinGene
##### BMI PRS
BMIprs_TG <- read.table("P:/IGEMS/IGEMS_Research/PRS_work/Data/BMI_Yengo2018_SBayesR_STR_TwinGene.txt", header=T, sep="\t")
##### Standardize and add study variable
BMIprs_TG$zPRS_BMI <- scale(BMIprs_TG$BMI_Yengo2018)
BMIprs_TG$study <- 'TwinGene'

########## SALT-Y
##### BMI PRS
BMIprs_SY <- read.table("P:/IGEMS/IGEMS_Research/PRS_work/Data/BMI_Yengo2018_SBayesR_STR_SALTY.txt", header=T, sep="\t")
##### Standardize and add study variable
BMIprs_SY$zPRS_BMI <- scale(BMIprs_SY$BMI_Yengo2018)
BMIprs_SY$study <- 'SALT-Y'

########## Stack together, and merge with the PCs
prs <- rbind(BMIprs_GOSH, BMIprs_TG, BMIprs_SY)
prs_pc <- merge(prs, pc, by=c("twinnr", "study"))
# n=19361

########## Import and add outcome data + PC adjusted PGS
BMI <- read_sas("P:/Dementia_IK/SfoEpi/GPdisc_CVD/Data/GOSH/gpdisc_cvd_str_adata_20221212.sas7bdat")
colnames(BMI)
BMI <- BMI[c(1,2,4,19,20,28)]
# n=17988

prs_out <- merge(prs_pc, BMI, by="twinnr")
# n=17988

########## Plot raw PRS values
boxplot(prs_out$BMI_Yengo2018 ~ prs_out$study, par(cex.axis=2.5), ylab="", xlab="BMI, Yengo 2018", par(cex.lab=3))

# Very nice and evenly distributed!
# Saving plot
pdf(file="P:/Dementia_IK/SfoEpi/GPdisc_CVD/Output/GOSH/PRS_BMI_posthocQC_20230123_box.pdf",width=15,height=10)
par(mfrow=c(1,1)) 
# Boxplots of the raw PGSs by study
boxplot(prs_out$BMI_Yengo2018 ~ prs_out$study, par(cex.axis=2.5), ylab="", xlab="BMI, Yengo 2018", par(cex.lab=3))
dev.off()

########## Compare the predictive ability of the PC adjusted PGS to that of the raw PGS adjusted for PCs
# Standardize the two scores
prs_out$z_PRS <- scale(prs_out$BMI_Yengo2018)
prs_out$z_PRSadj <- scale(prs_out$PRS_BMI_PCadj)

# 1. Phe standardized PGS, additionally adjusted for PCs
lr_raw0 <- lm(BMI~ sex + BMI_age + study + PC1 + PC2 + PC3 + PC4 + PC5, data=prs_out)
lr_raw1 <-lm(BMI~ z_PRS + sex + BMI_age + study + PC1 + PC2 + PC3 + PC4 + PC5, data=prs_out)
# Calculate variance explained 
summary(lr_raw1)$r.squared - summary(lr_raw0)$r.squared
# R2 = 0.1086014
# Get coefficients
vcovCL <- cluster.vcov(lr_raw1, prs_out$pairid) 
coeftest(lr_raw1, vcovCL)
#                 Estimate Std. Error  t value  Pr(>|t|)    
# (Intercept)   23.1962434  0.1918927 120.8813 < 2.2e-16 ***
# z_PRS          1.1165031  0.0265537  42.0470 < 2.2e-16 ***
# sex           -0.9800408  0.0493595 -19.8551 < 2.2e-16 ***
# BMI_age        0.0616521  0.0031524  19.5573 < 2.2e-16 ***
# studySALT-Y   -0.0083009  0.0909270  -0.0913 0.9272615    
# studyTwinGene -0.0776330  0.0871006  -0.8913 0.3727787    
# PC1            9.1159944  2.3544513   3.8718 0.0001084 ***
# PC2           -6.4629514  2.4624747  -2.6246 0.0086831 ** 
# PC3           -4.2643601  2.7228672  -1.5661 0.1173362    
# PC4            1.4029122  2.6618244   0.5270 0.5981660    
# PC5            0.0606915  2.8130481   0.0216 0.9827872

# 2. The PC adjusted (and standardized) PGS
lr_adj0 <- lm(BMI~ sex + BMI_age + study, data=prs_out)
lr_adj1 <-lm(BMI~ z_PRSadj + sex + BMI_age + study, data=prs_out)
# Calculate variance explained 
summary(lr_adj1)$r.squared - summary(lr_adj0)$r.squared
# R2 = 0.1084947
# Get coefficients
vcovCL <- cluster.vcov(lr_adj1, prs_out$pairid) 
coeftest(lr_adj1, vcovCL)
#                 Estimate Std. Error  t value  Pr(>|t|)    
# (Intercept)   23.231640   0.192068 120.9555   <2e-16 ***
# z_PRSadj       1.110139   0.026465  41.9474   <2e-16 ***
# sex           -0.979271   0.049499 -19.7836   <2e-16 ***
# BMI_age        0.061104   0.003163  19.3181   <2e-16 ***
# studySALT-Y    0.003143   0.091056   0.0345   0.9725    
# studyTwinGene -0.078180   0.087198  -0.8966   0.3700


########## Done!
