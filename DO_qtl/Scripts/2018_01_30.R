#############################################################
# Mediation Analysis for Ins Sec QTL Peaks
# January 30th, 2018
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

# table with old and new phenotype names
renaming <- read.csv("QTL_parameters/18_01_24_ins_renaming_parameters.csv", stringsAsFactors=FALSE)
stopifnot( all(!is.na(match(renaming$old_name, colnames(pheno)))) )

# re-name the phenotypes
phe_col <- colnames(pheno)
phe_col[match(renaming$old_name, phe_col)] <- renaming$new_name
colnames(pheno) <- phe_col

# Read Ins. Sec. QTL peak summary CSV
QTLsummary = read.csv("Derived_data/2018_01_29_ins_sec_qtl_summary_thresh_6.csv", stringsAsFactors = FALSE)

# Functions for mediation analysis
sapply(list.files("Scripts/Mediation_Tools", full.names = TRUE), source)
load("Data/mouse.chrlen.rda")

# Set mediator to gene expression datafile 
mediator = rankz.mrna
# Set annotation to gene expression annotation file
annotation = annot.mrna
names(annotation)[colnames(annotation) == "middle_point"] = "pos"
# Set covar to covariate model matrix
pheno$sex = as.factor(pheno$sex)
pheno$DOwave = as.factor(pheno$DOwave)
covar = model.matrix(~sex + diet_days + DOwave, 
                     pheno[rownames(mediator),])[,-1,drop=FALSE]

stopifnot(all(sub("_log.*", "", QTLsummary$analyte) %in% colnames(pheno)))
out = vector("list", nrow(QTLsummary))
names(out) = QTLsummary$analyte
for(i in 1:length(out))
{
  cat("Batch", i, "of", length(out), "\n")
  target = sub("_log.*", "", QTLsummary$analyte[i])
  target = log(pheno[rownames(mediator),
                     target, 
                     drop = FALSE])
  
  chr = QTLsummary$chr[i]
  pos = QTLsummary$pos[i]
  chr_pos = paste(chr, pos*1e+06, sep = "_")
  
  qtl.geno = genoprobs[[chr]][rownames(mediator),,chr_pos]
  
  med <- mediation.scan(target = target,
                        mediator = mediator,
                        annotation = annotation,
                        covar = covar,
                        qtl.geno = qtl.geno)
  med$chr[med$chr == "MT"] = "M"
  
  out[[i]] = med
}

for(i in 31:63)
{
  quartz()
  plot(out[[i]])
}

x =scan1(genoprobs, pheno = target, K, addcovar=covar, cores=0)
par(mfrow = c(2,1))
plot(x,lodcolumn = 1, map = map)
plot(y,lodcolumn = 1, map = map)


