#############################################################
library(qtl2)
library(qtl2convert)
library(tidyverse)
setwd("/Users/s-ishimt/Desktop/Academic_Internship/DO_qtl")
#############################################################
# load data file
load("Data/Attie_islet_secr_data_v1.Rdata")
ins_sec_qtl = readRDS("Derived_data/2018_01_02_ins_output.rds")
ins_sec_peaks = read.csv("Derived_data/2018_01_05_ins_sec_qtl_summary_thresh_7.csv")
renaming <- read.csv("QTL_parameters/18_01_02_ins_renaming_parameters.csv", stringsAsFactors=FALSE)

# filter peaks summary file to contain target chr 11 peaks
ins_sec_peaks = ins_sec_peaks[order(ins_sec_peaks$chr), ]
ins_sec_chr1_peaks = ins_sec_peaks %>%
  filter(chr == 1)
ins_sec_chr1_peaks$analyte = c("G33_log_covlogIns_per_islet", "GLP1_log_covlogIns_per_islet","AA_log_covlogIns_per_islet")

# create parameters file for effect plots
output = as.character(ins_sec_chr1_peaks$analyte)
phenotypes = NA
phenotypes[grep("_log", output)] = gsub("_log.*", "", output[grep("_log", output)])
phenotypes[grep("_rankz", output)] = gsub("_rankz.*", "", output[grep("_rankz", output)])
transformation = NA
transformation[grep("_log", output)] = "log"
transformation[grep("_rankz", output)] = "qtl::nqrank"
covariates = rep("sex;diet_days;DOwave", length(output))
covariates[grep("_cov", output)] = paste(
  "sex;diet_days;DOwave", gsub(".*_covlog","log_", output[grep("_cov", output)]), sep = ";")

analyses = cbind(output, phenotypes, transformation, covariates, ins_sec_chr1_peaks[,c("chr", "pos")], stringsAsFactors = FALSE)

# re-name the phenotypes
phe_col <- colnames(pheno)
phe_col[match(renaming$old_name, phe_col)] <- renaming$new_name
colnames(pheno) <- phe_col

# covariate data frame
# grab names
covar_names <- unique(unlist(strsplit(analyses$covariates, ";")))
covar_transf = sub("log_", "", covar_names)
covar_transf = c("sex", "diet_days", "DOwave", "Ins_per_islet")
# pull out covariate data, adding columns (Glu0 - 50) + DOwave3,4,5
covar <- pheno[,covar_transf]

# take logs of appropriate columns
covar[,grep("^log", covar_names)] <- log(covar[,grep("^log", covar_names)]) 
# add "log_" to front of those column names
colnames(covar)[grep("^log", covar_names)] <- paste0("log_", colnames(covar)[grep("^log", covar_names)])

# do genome scans with row of analysese
out <- vector("list", nrow(analyses))
names(out) <- analyses$output

for(i in 1:nrow(analyses)){
  cat("Batch", i, "of", nrow(analyses), "\n")
  
  # pull out the target phenotype
  this_pheno = pheno[,analyses$pheno[i],drop=FALSE]
  
  # apply the transformations
  f_name <- analyses$transformation[i]
  if(f_name == "qtl::nqrank")
  {
    f = qtl::nqrank
  } else
  {
    f = get(f_name)
  }
  this_pheno = f(this_pheno)
  
  # paste in output name as phenotype name
  colnames(this_pheno) <- analyses$output[i]
  
  # construct covariate matrix (omitting the intercept)
  covarX <- model.matrix(as.formula(
    paste0("~", gsub(";", "+", analyses$covariates[i]))), covar)[,-1,drop=FALSE]
  
  # set chromosome
  chr = analyses$chr[i]
  
  # run scan1blup
  out[[i]] <- scan1blup(genoprobs[,chr], this_pheno, K[[chr]], addcovar=covarX, cores=0)
}

# save output to RDS file
saveRDS(out, "Derived_data/2018_01_16/scan1blups_chr1.rds")

# plot and save in pdf
names(out) = c("G33_log_covlogInsPerislet", "GLP1_log_covlogInsPerIslet", "AA_log_covlogInsPerIslet")
pdf(file = "Plots/2018_01_16/chr1_qtlblups.pdf", width = 10, height = 6)
for(i in seq_along(out)){
  chr = as.numeric(sub("_.*", "", rownames(out[[i]])[1]))
  if(i %in% grep("_cov",names(out))){
    cov_name = unlist(strsplit(names(out)[i],"_cov"))[2]
    plot_name = unlist(strsplit(names(out)[i],"_cov"))[1]
    plot(out[[i]], map = map[[chr]], columns = 1:8, col = CCcolors, 
         main = plot_name, scan1_output = ins_sec_qtl[,names(out)[i], drop = FALSE])
    mtext(paste("Covariate:", cov_name), cex = 1, line = 0.35)
  }else{
    plot(out[[i]], map = map[[chr]], columns = 1:8, col = CCcolors,
         main = names(out)[i], scan1_output = ins_sec_qtl[,names(out)[i], drop = FALSE])
  }
}
dev.off()


