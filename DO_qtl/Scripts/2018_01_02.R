#############################################################
# Analysis of insulin secretion phenotypes for 2nd QTL paper
# January 2nd, 2018
# - Taki Ishimura
#############################################################
library(qtl2)
setwd("/Users/s-ishimt/Desktop/Academic_Internship/DO_qtl")
#############################################################

# load data file
load("Data/Attie_islet_secr_data_v1.Rdata")
ls()
# contains:
# genome.build  - GRCm38
# genoprobs     - R/qtl2 genoprobs on 69k grid
# K             - list of "loco" kinship matrices (500 x 500)
# map           - physical map for the 69k grid
# markers       - data frame for 69k grid with marker, chr, pos, cM, bp
# pheno         - phenotypes (483 x 57)
# pheno_dict    - phenotype dictionary
# annot.mrna    - 21,771 x 9 mRNA annotations (id,symbol,chr,start,end,strand,middle_point,nearest_marker,biotype)
# expr.mrna     - 378 x 21,771 mRNA phenotypes
# rankz.mrna    - 378 x 21,771 mRNA phenotypes (rankZ)
# raw.mrna      - 378 x 21,771 mRNA phenotypes (raw)

# table with old and new phenotype names
renaming <- read.csv("QTL_parameters/18_01_02_ins_renaming_parameters.csv", stringsAsFactors=FALSE)
stopifnot( all(!is.na(match(renaming$old_name, colnames(pheno)))) )

# re-name the phenotypes
phe_col <- colnames(pheno)
phe_col[match(renaming$old_name, phe_col)] <- renaming$new_name
colnames(pheno) <- phe_col

# table with analyses
analyses <- read.csv("QTL_parameters/18_01_02_ins_analyses_parameters.csv", stringsAsFactors=FALSE)
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
# add "log_" to front of those column names
colnames(covar)[grep("^log", covar_names)] <- paste0("log_", colnames(covar)[grep("^log", covar_names)])

# unique sets of covariates
u_covar <- unique(analyses$covar)

# do genome scans with each unique set of covariates
out <- vector("list", length(u_covar))
names(out) <- u_covar
for(i in seq_along(u_covar)) {
  cat("Batch", i, "of", length(u_covar), "\n")
  # pull out analysis rows with this set of covariates
  these_rows <- which(analyses$covar == u_covar[i])
  
  # pull out this set of phenotypes
  these_pheno <- pheno[,analyses$pheno[these_rows],drop=FALSE]
  
  # apply the transformations
  for(j in 1:ncol(these_pheno)) {
    f_name <- analyses$transformation[these_rows[j]]
    if(f_name=="qtl::nqrank") f <- qtl::nqrank
    else f <- get(f_name)
    these_pheno[,j] <- f(these_pheno[,j])
  }
  
  # paste in output name as phenotype name
  colnames(these_pheno) <- analyses$output[these_rows]
  
  # construct covariate matrix (omitting the intercept)
  covarX <- model.matrix(as.formula(paste0("~", gsub(";", "+", u_covar[i]))), covar)[,-1,drop=FALSE]
  
  out[[i]] <- scan1(genoprobs, these_pheno, K, addcovar=covarX, cores=0)
}

# combine output into a single scan1 object with cbind()
out <- do.call("cbind", out)

# reorder the columns in the results (oog...need a helper function for this)
out_rev <- out[,analyses$output]
attr(out_rev, "sample_size") <- attr(out, "sample_size")[analyses$output]
attr(out_rev, "hsq") <- attr(out, "hsq")[,analyses$output]
class(out_rev) <- class(out)
out <- out_rev

# save output to RDS file
saveRDS(out, "Derived_data/2018_01_02_ins_output.rds")


# print QTL plots
pdf(file = "Plots/2018_01_02/ins_sec_qtl.pdf", width = 8, height = 4)
for(i in 1:ncol(out)){
  if(i %in% grep("_cov",colnames(out))){
    cov_name = unlist(strsplit(colnames(out)[i],"_cov"))[2]
    plot_name = unlist(strsplit(colnames(out)[i],"_cov"))[1]
    plot(out,lodcolumn = i, map = map, main = plot_name)
    mtext(paste("Covariate:", cov_name), cex = 1, line = 0.35)
  }else{
    plot(out,lodcolumn = i, map = map, main = colnames(out)[i])
  }
}
dev.off()
