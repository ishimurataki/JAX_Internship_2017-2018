#############################################################
# Interactive Sex Effects
# February 5th, 2018
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
pheno$sex = as.factor(pheno$sex)
pheno$DOwave = as.factor(pheno$DOwave)

# table with old and new phenotype names
renaming <- read.csv("QTL_parameters/18_01_24_ins_renaming_parameters.csv", stringsAsFactors=FALSE)
stopifnot( all(!is.na(match(renaming$old_name, colnames(pheno)))) )

# re-name the phenotypes
phe_col <- colnames(pheno)
phe_col[match(renaming$old_name, phe_col)] <- renaming$new_name
colnames(pheno) <- phe_col

# table with analyses
analyses <- read.csv("QTL_parameters/18_01_24_ins_analyses_parameters.csv", stringsAsFactors=FALSE)
stopifnot( all(!is.na(match(analyses$pheno, colnames(pheno)))) )

# covariate data frame
# grab names
covar_names <- unique(unlist(strsplit(analyses$covar, ";")))
# strip off transformation
covar_notransf <- sub("^log_", "", covar_names)
# pull out covariate data, adding columns (Glu0 - 50) + DOwave3,4,5
covar <- cbind(pheno,
               Glu0m50=pheno$Glu0 - 50,
               DOwave3=(pheno$DOwave==3)*1,
               DOwave4=(pheno$DOwave==4)*1,
               DOwave5=(pheno$DOwave==5)*1)[,covar_notransf]
# take logs of appropriate columns
covar[,grep("^log", covar_names)] <- log(covar[,grep("^log", covar_names)]) ##### NOT ACTUALLY TAKING THE LOG OF THE COVARIATES ##########
covar$sex = as.factor(covar$sex)
covar$DOwave = as.factor(covar$DOwave)
# add "log_" to front of those column names
colnames(covar)[grep("^log", covar_names)] <- paste0("log_", colnames(covar)[grep("^log", covar_names)])

# unique sets of covariates
u_covar <- unique(analyses$covar)

# do genome scans with each unique set of covariates
out <- vector("list", 2*nrow(analyses))
outnames = vector("character")
for(i in 1:nrow(analyses)){
  add = paste("add", analyses$output[i], sep = "_")
  int = paste("int", analyses$output[i], sep = "_")
  add_int = c(add, int)
  outnames = c(outnames, add_int)
}
names(out) = outnames

for(i in 1:(length(out)/2)){
  cat("Batch", i, "of", (length(out)/2), "\n")
  # pull out analysis rows with this set of covariates
  phenotype = pheno[,analyses$pheno[i], drop = FALSE]
  
  f_name = analyses$transformation[i]
  f = get(f_name)
  phenotype = f(phenotype)
  
  # construct covariate matrix (omitting the intercept)
  covarX = model.matrix(as.formula(paste0("~", gsub(";", "+", analyses$covariates[i]))), 
                         covar)[,-1,drop=FALSE]
  covarSex = model.matrix(~sex, covar)[,-1,drop=FALSE]
  
  additive = scan1(genoprobs, phenotype, K, addcovar=covarX, cores=0)
  
  interactive = scan1(genoprobs, phenotype, K, addcovar=covarX, intcovar = covarSex, cores=0)
  
  colnames(additive) = paste("additive", names(phenotype), sep = "_")
  colnames(interactive) = paste("interactive", names(phenotype), sep = "_")
  
  out[[2*i-1]] = additive
  out[[2*i]] = interactive
}

pdf(file = "Plots/2018_02_06/sex_effects.pdf", width = 8, height = 7.5)
for(i in 1:(ncol(testout)/2)){
  print(i)
  par(mfrow = c(3,1), 
      oma = c(5,4,4,4) + 0.1,
      mar = c(1.8,1,1.5,1) + 0.1)
  
  additive = out[[2*i-1]]
  interactive = out[[2*i]]
  difference = out[[2*i]] - out[[2*i-1]]
  
  plot(additive, map = map, main = "sex additive")
  plot(interactive, map = map, main = "sex interactive")
  plot(difference, map = map, main = "interactive - additive")
  title(analyses$output[i], outer = TRUE, cex.main = 1.7)
}
dev.off()
