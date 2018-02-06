#############################################################
# Sex-interactive QTL for Gcg_secreted cov = Gcg_content
# January 26th, 2018
# - Taki Ishimura
#############################################################
library(ggplot2)
setwd("/Users/s-ishimt/Desktop/Academic_Internship/DO_qtl")
#############################################################

# load data file
load("Data/Attie_islet_secr_data_v1.Rdata")
load("Data/pheno_clin_v6.RData")
pheno = pheno_clin
pheno_dict = pheno_clin_dict
rm(pheno_clin, pheno_clin_dict)

Gcg_secreted = log(pheno[,"Gcg_secreted", drop = FALSE])
covar = model.matrix(~sex + DOwave + diet_days, data=pheno)[,-1]
intcovar = model.matrix(~sex, data = pheno)[,-1]


additive = scan1(genoprobs, pheno = Gcg_secreted, K, addcovar=covar, cores=0)
interactive = scan1(genoprobs, pheno = Gcg_secreted, K, addcovar=covar, cores=0, intcovar = intcovar)

pdf(file = "Plots/2018_01_26/sex_interactive_additive.pdf", width = 10, height = 10)
par(mfrow = c(3,1))
plot(additive,lodcolumn = 1, map = map, 
     main = "Gcg_secreted|Gcg_content, sex = additive")
plot(interactive,lodcolumn = 1, map = map, 
     main = "Gcg_secreted|Gcg_content, sex = interactive")
plot(interactive - additive,lodcolumn = 1, map = map, 
     main = "interactive - additive")
dev.off()
