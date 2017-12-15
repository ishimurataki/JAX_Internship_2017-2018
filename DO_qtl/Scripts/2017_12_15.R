###########################################
#  Attie DO Islet Secretion Data Analysis
#  QC and scratch pad
#  GAC Decemeber 2017
# -Taki Ishimura
###########################################
library(tidyverse)
library(qtl2)
library(qtl2convert)

setwd("/Users/s-ishimt/Desktop/Academic_Internship/DO_qtl")
###########################################
# Set Up

####
# load data
load("Data/Attie_islet_secr_data_v1.Rdata")
ls()
# "annot.mrna"   "expr.mrna"    "genome.build" "genoprobs"    "K"
# "map"          "markers"      "pheno"        "pheno_dict"   "rankz.mrna"  "raw.mrna"

####
# convert pheno to a tibble save as Pheno
Pheno <- as_data_frame(pheno) %>%
  select(mouse, birthdate, DOwave, sex, chrM, chrY, diet_days, food_ave,
         num_islets, Ins_per_islet, WPIC, Glu_0min, Ins_0min, HOMA_B_0min,
         HOMA_IR_0min, Glu_tAUC, Glu_iAUC, Ins_tAUC, Ins_iAUC, Glu_sac,
         Ins_sac, TG_sac, weight_sac, Gcg_content, Gcg_secreted, G33_ins_secrete_gm,
         G83_ins_secrete_gm, G167_ins_secrete_gm, KCl_G33_ins_secrete_gm,
         GLP1_G83_ins_secrete_gm, AA_G83_ins_secrete_gm, PA_G167_ins_secrete_gm ) %>%
  rename(Glu0 = Glu_0min,
         Ins0 = Ins_0min,
         HOMA_B = HOMA_B_0min,
         HOMA_IR = HOMA_IR_0min,
         Glu = Glu_sac, Ins = Ins_sac,
         TG = TG_sac, weight = weight_sac,
         G33 = G33_ins_secrete_gm,
         G83 = G83_ins_secrete_gm,
         G167 = G167_ins_secrete_gm,
         KCl = KCl_G33_ins_secrete_gm,
         GLP1 = GLP1_G83_ins_secrete_gm,
         AA = AA_G83_ins_secrete_gm,
         PA = PA_G167_ins_secrete_gm)
names(Pheno)
# "mouse"         "birthdate"     "DOwave"        "sex"           "chrM"
# "chrY"          "diet_days"     "food_ave"      "num_islets"    "Ins_per_islet"
# "WPIC"          "Glu0"          "Ins0"          "HOMA_B"        "HOMA_IR"
# "Glu_tAUC"      "Glu_iAUC"      "Ins_tAUC"      "Ins_iAUC"      "Glu"
# "Ins"           "TG"            "weight"        "Gcg_content"   "Gcg_secreted"
# "G33"           "G83"           "G167"          "KCl"           "GLP1"
# "AA"            "PA"

#time to get scanning. 
plot_location = "/Users/s-ishimt/Desktop/Academic_Internship/DO_qtl/Plots/2017_12_15/"

# standard baseline covariates are sex, wave, and diet_days.
# Phenotype: num_islets, Transformation: sqrt transform, Additional Covariates: None.
phenotype = "num_islets"
transformation = sqrt
pdfname = "num_islets_sqrt"

plotname = paste(phenotype, "sqrt", sep = "_")
phenotype = transformation(as.data.frame(Pheno)[,phenotype,drop = FALSE])
addcovar = model.matrix(~sex + DOwave + diet_days, data = Pheno)

qtl = scan1(genoprobs = genoprobs, pheno = phenotype,
            kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)

filename = paste(plot_location, pdfname, ".pdf", sep = "")
pdf(file = filename, width = 8, height = 4)
plot(qtl, map = map, main = plotname)
dev.off()

# Phenotype: Ins_per_islet, Transformation: log, Additional Covariates: None.
phenotype = "Ins_per_islet"
transformation = log
pdfname = "Ins_per_islet_log"

plotname = paste(phenotype, "log", sep = "_")
phenotype = transformation(as.data.frame(Pheno)[,phenotype,drop = FALSE])
addcovar = model.matrix(~sex + DOwave + diet_days, data = Pheno)

qtl = scan1(genoprobs = genoprobs, pheno = phenotype,
            kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)

filename = paste(plot_location, pdfname, ".pdf", sep = "")
pdf(file = filename, width = 8, height = 4)
plot(qtl, map = map, main = plotname)
dev.off()

# Phenotype: WPIC, Transformation: log, Additional Covariates: None.
phenotype = "WPIC"
transformation = log
pdfname = "WPIC_log"

plotname = paste(phenotype, "log", sep = "_")
phenotype = transformation(as.data.frame(Pheno)[,phenotype,drop = FALSE])
addcovar = model.matrix(~sex + DOwave + diet_days, data = Pheno)

qtl = scan1(genoprobs = genoprobs, pheno = phenotype,
            kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)

filename = paste(plot_location, pdfname, ".pdf", sep = "")
pdf(file = filename, width = 8, height = 4)
plot(qtl, map = map, main = plotname)
dev.off()


# Phenotype: Glu0, Transformation: log, Additional Covariates: None.
phenotype = "Glu0"
transformation = log
pdfname = "Glu0_log"

plotname = paste(phenotype, "log", sep = "_")
phenotype = transformation(as.data.frame(Pheno)[,phenotype,drop = FALSE])
addcovar = model.matrix(~sex + DOwave + diet_days, data = Pheno)

qtl = scan1(genoprobs = genoprobs, pheno = phenotype,
            kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)

filename = paste(plot_location, pdfname, ".pdf", sep = "")
pdf(file = filename, width = 8, height = 4)
plot(qtl, map = map, main = plotname)
dev.off()

# Phenotype: Glu0, Transformation: log, Additional Covariates: Ins0
phenotype = "Glu0"
transformation = log
pdfname = "Glu0_log_COVIns0"

plotname = paste(phenotype, "log", sep = "_")
phenotype = transformation(as.data.frame(Pheno)[,phenotype,drop = FALSE])
addcovar = model.matrix(~sex + DOwave + diet_days + log(Ins0), data = Pheno)

qtl = scan1(genoprobs = genoprobs, pheno = phenotype,
            kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)

filename = paste(plot_location, pdfname, ".pdf", sep = "")
pdf(file = filename, width = 8, height = 4)
plot(qtl, map = map, main = plotname)
mtext(paste("Covariate:", "Ins0_log"), cex = 1, line = 0.35)
dev.off()

# Phenotype: Ins0, Transformation: log, Additional Covariates: None
phenotype = "Ins0"
transformation = log
pdfname = "Ins0_log"

plotname = paste(phenotype, "log", sep = "_")
phenotype = transformation(as.data.frame(Pheno)[,phenotype,drop = FALSE])
addcovar = model.matrix(~sex + DOwave + diet_days, data = Pheno)

qtl = scan1(genoprobs = genoprobs, pheno = phenotype,
            kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)

filename = paste(plot_location, pdfname, ".pdf", sep = "")
pdf(file = filename, width = 8, height = 4)
plot(qtl, map = map, main = plotname)
dev.off()

# Phenotype: HOMA_B, Transformation: log, Additional Covariates: None
phenotype = "HOMA_B"
transformation = log
pdfname = "HOMA_B_log"

plotname = paste(phenotype, "log", sep = "_")
phenotype = transformation(as.data.frame(Pheno)[,phenotype,drop = FALSE])
addcovar = model.matrix(~sex + DOwave + diet_days, data = Pheno)

qtl = scan1(genoprobs = genoprobs, pheno = phenotype,
            kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)

filename = paste(plot_location, pdfname, ".pdf", sep = "")
pdf(file = filename, width = 8, height = 4)
plot(qtl, map = map, main = plotname)
dev.off()

# Phenotype: HOMA_IR, Transformation: log, Additional Covariates: None
phenotype = "HOMA_IR"
transformation = log
pdfname = "HOMA_IR_log"

plotname = paste(phenotype, "log", sep = "_")
phenotype = transformation(as.data.frame(Pheno)[,phenotype,drop = FALSE])
addcovar = model.matrix(~sex + DOwave + diet_days, data = Pheno)

qtl = scan1(genoprobs = genoprobs, pheno = phenotype,
            kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)

filename = paste(plot_location, pdfname, ".pdf", sep = "")
pdf(file = filename, width = 8, height = 4)
plot(qtl, map = map, main = plotname)
dev.off()

# Phenotype: Glu_tAUC, Transformation: log, Additional Covariates: None
phenotype = "Glu_tAUC"
transformation = log
pdfname = "Glu_tAUC_log"

plotname = paste(phenotype, "log", sep = "_")
phenotype = transformation(as.data.frame(Pheno)[,phenotype,drop = FALSE])
addcovar = model.matrix(~sex + DOwave + diet_days, data = Pheno)

qtl = scan1(genoprobs = genoprobs, pheno = phenotype,
            kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)

filename = paste(plot_location, pdfname, ".pdf", sep = "")
pdf(file = filename, width = 8, height = 4)
plot(qtl, map = map, main = plotname)
dev.off()

# Phenotype: Glu_tAUC, Transformation: log, Additional Covariates: Glu0
phenotype = "Glu_tAUC"
transformation = log
pdfname = "Glu_tAUC_log_COVGlu0_log"

plotname = paste(phenotype, "log", sep = "_")
phenotype = transformation(as.data.frame(Pheno)[,phenotype,drop = FALSE])
addcovar = model.matrix(~sex + DOwave + diet_days + log(Glu0), data = Pheno)

qtl = scan1(genoprobs = genoprobs, pheno = phenotype,
            kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)

filename = paste(plot_location, pdfname, ".pdf", sep = "")
pdf(file = filename, width = 8, height = 4)
plot(qtl, map = map, main = plotname)
mtext(paste("Covariate:", "Glu0_log"), cex = 1, line = 0.35)
dev.off()

# Phenotype: = Glu_iAUC, Transformation: log, Additional Covariates: None
phenotype = "Glu_iAUC"
transformation = log
pdfname = "Glu_iAUC_log"

plotname = paste(phenotype, "log", sep = "_")
phenotype = transformation(as.data.frame(Pheno)[,phenotype,drop = FALSE])
addcovar = model.matrix(~sex + DOwave + diet_days, data = Pheno)

qtl = scan1(genoprobs = genoprobs, pheno = phenotype,
            kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)

filename = paste(plot_location, pdfname, ".pdf", sep = "")
pdf(file = filename, width = 8, height = 4)
plot(qtl, map = map, main = plotname)
dev.off()

# Phenotype: = Ins_tAUC, Transformation: log, Additional Covariates: None
phenotype = "Ins_tAUC"
transformation = log
pdfname = "Ins_tAUC_log"

plotname = paste(phenotype, "log", sep = "_")
phenotype = transformation(as.data.frame(Pheno)[,phenotype,drop = FALSE])
addcovar = model.matrix(~sex + DOwave + diet_days, data = Pheno)

qtl = scan1(genoprobs = genoprobs, pheno = phenotype,
            kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)

filename = paste(plot_location, pdfname, ".pdf", sep = "")
pdf(file = filename, width = 8, height = 4)
plot(qtl, map = map, main = plotname)
dev.off()

# Phenotype: = Ins_tAUC, Transformation: log, Additional Covariates: Ins0
phenotype = "Ins_tAUC"
transformation = log
pdfname = "Ins_tAUC_log_COVIns0_log"

plotname = paste(phenotype, "log", sep = "_")
phenotype = transformation(as.data.frame(Pheno)[,phenotype,drop = FALSE])
addcovar = model.matrix(~sex + DOwave + diet_days + log(Ins0), data = Pheno)

qtl = scan1(genoprobs = genoprobs, pheno = phenotype,
            kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)

filename = paste(plot_location, pdfname, ".pdf", sep = "")
pdf(file = filename, width = 8, height = 4)
plot(qtl, map = map, main = plotname)
mtext(paste("Covariate:", "Ins0_log"), cex = 1, line = 0.35)
dev.off()

# Phenotype: = Ins_iAUC, Transformation: log, Additional Covariates: None
phenotype = "Ins_iAUC"
transformation = log
pdfname = "Ins_iAUC_log"

plotname = paste(phenotype, "log", sep = "_")
phenotype = transformation(as.data.frame(Pheno)[,phenotype,drop = FALSE])
addcovar = model.matrix(~sex + DOwave + diet_days, data = Pheno)

qtl = scan1(genoprobs = genoprobs, pheno = phenotype,
            kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)

filename = paste(plot_location, pdfname, ".pdf", sep = "")
pdf(file = filename, width = 8, height = 4)
plot(qtl, map = map, main = plotname)
dev.off()

# Phenotype: = Glu, Transformation: log, Additional Covariates: None
phenotype = "Glu"
transformation = log
pdfname = "Glu_log"

plotname = paste(phenotype, "log", sep = "_")
phenotype = transformation(as.data.frame(Pheno)[,phenotype,drop = FALSE])
addcovar = model.matrix(~sex + DOwave + diet_days, data = Pheno)

qtl = scan1(genoprobs = genoprobs, pheno = phenotype,
            kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)

filename = paste(plot_location, pdfname, ".pdf", sep = "")
pdf(file = filename, width = 8, height = 4)
plot(qtl, map = map, main = plotname)
dev.off()

# Phenotype: = Ins, Transformation: log, Additional Covariates: None
phenotype = "Ins"
transformation = log
pdfname = "Ins_log"

plotname = paste(phenotype, "log", sep = "_")
phenotype = transformation(as.data.frame(Pheno)[,phenotype,drop = FALSE])
addcovar = model.matrix(~sex + DOwave + diet_days, data = Pheno)

qtl = scan1(genoprobs = genoprobs, pheno = phenotype,
            kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)

filename = paste(plot_location, pdfname, ".pdf", sep = "")
pdf(file = filename, width = 8, height = 4)
plot(qtl, map = map, main = plotname)
dev.off()

# Phenotype: = TG, Transformation: log, Additional Covariates: None
phenotype = "TG"
transformation = log
pdfname = "TG_log"

plotname = paste(phenotype, "log", sep = "_")
phenotype = transformation(as.data.frame(Pheno)[,phenotype,drop = FALSE])
addcovar = model.matrix(~sex + DOwave + diet_days, data = Pheno)

qtl = scan1(genoprobs = genoprobs, pheno = phenotype,
            kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)

filename = paste(plot_location, pdfname, ".pdf", sep = "")
pdf(file = filename, width = 8, height = 4)
plot(qtl, map = map, main = plotname)
dev.off()

# Phenotype: = Gcg_content, Transformation: log, Additional Covariates: None
phenotype = "Gcg_content"
transformation = log
pdfname = "Gcg_content_log"

plotname = paste(phenotype, "log", sep = "_")
phenotype = transformation(as.data.frame(Pheno)[,phenotype,drop = FALSE])
addcovar = model.matrix(~sex + DOwave + diet_days, data = Pheno)

qtl = scan1(genoprobs = genoprobs, pheno = phenotype,
            kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)

filename = paste(plot_location, pdfname, ".pdf", sep = "")
pdf(file = filename, width = 8, height = 4)
plot(qtl, map = map, main = plotname)
dev.off()

# Phenotype: = Gcg_secreted, Transformation: log, Additional Covariates: None
phenotype = "Gcg_secreted"
transformation = log
pdfname = "Gcg_secreted_log"

plotname = paste(phenotype, "log", sep = "_")
phenotype = transformation(as.data.frame(Pheno)[,phenotype,drop = FALSE])
addcovar = model.matrix(~sex + DOwave + diet_days, data = Pheno)

qtl = scan1(genoprobs = genoprobs, pheno = phenotype,
            kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)

filename = paste(plot_location, pdfname, ".pdf", sep = "")
pdf(file = filename, width = 8, height = 4)
plot(qtl, map = map, main = plotname)
dev.off()

# Phenotype: = Gcg_secreted, Transformation: log, Additional Covariates: Gcg_content
phenotype = "Gcg_secreted"
transformation = log
pdfname = "Gcg_secreted_log_COVGcg_content_log"

plotname = paste(phenotype, "log", sep = "_")
phenotype = transformation(as.data.frame(Pheno)[,phenotype,drop = FALSE])
addcovar = model.matrix(~sex + DOwave + diet_days + log(Gcg_content), data = Pheno)

qtl = scan1(genoprobs = genoprobs, pheno = phenotype,
            kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)

filename = paste(plot_location, pdfname, ".pdf", sep = "")
pdf(file = filename, width = 8, height = 4)
plot(qtl, map = map, main = plotname)
mtext(paste("Covariate:", "Gcg_content_log"), cex = 1, line = 0.35)
dev.off()

# Phenotype: = G33, Transformation: log, Additional Covariates: None
phenotype = "G33"
transformation = log
pdfname = "G33_log"

plotname = paste(phenotype, "log", sep = "_")
phenotype = transformation(as.data.frame(Pheno)[,phenotype,drop = FALSE])
addcovar = model.matrix(~sex + DOwave + diet_days, data = Pheno)

qtl = scan1(genoprobs = genoprobs, pheno = phenotype,
            kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)

filename = paste(plot_location, pdfname, ".pdf", sep = "")
pdf(file = filename, width = 8, height = 4)
plot(qtl, map = map, main = plotname)
dev.off()

# Phenotype: G83, Transformation: log, Additional Covariates: None
phenotype = "G83"
transformation = log
pdfname = "G83_log"

plotname = paste(phenotype, "log", sep = "_")
phenotype = transformation(as.data.frame(Pheno)[,phenotype,drop = FALSE])
addcovar = model.matrix(~sex + DOwave + diet_days, data = Pheno)

qtl = scan1(genoprobs = genoprobs, pheno = phenotype,
            kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)

filename = paste(plot_location, pdfname, ".pdf", sep = "")
pdf(file = filename, width = 8, height = 4)
plot(qtl, map = map, main = plotname)
dev.off()


# Phenotype: G167, Transformation: log, Additional Covariates: None
phenotype = "G167"
transformation = log
pdfname = "G167_log"

plotname = paste(phenotype, "log", sep = "_")
phenotype = transformation(as.data.frame(Pheno)[,phenotype,drop = FALSE])
addcovar = model.matrix(~sex + DOwave + diet_days, data = Pheno)

qtl = scan1(genoprobs = genoprobs, pheno = phenotype,
            kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)

filename = paste(plot_location, pdfname, ".pdf", sep = "")
pdf(file = filename, width = 8, height = 4)
plot(qtl, map = map, main = plotname)
dev.off()

# Phenotype: G167, Transformation: log, Additional Covariates: G33
phenotype = "G167"
transformation = log
pdfname = "G167_log_COVG33_log"

plotname = paste(phenotype, "log", sep = "_")
phenotype = transformation(as.data.frame(Pheno)[,phenotype,drop = FALSE])
addcovar = model.matrix(~sex + DOwave + diet_days + log(G33), data = Pheno)

qtl = scan1(genoprobs = genoprobs, pheno = phenotype,
            kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)

filename = paste(plot_location, pdfname, ".pdf", sep = "")
pdf(file = filename, width = 8, height = 4)
plot(qtl, map = map, main = plotname)
mtext(paste("Covariate:", "G33_log"), cex = 1, line = 0.35)
dev.off()

# Phenotype: G167, Transformation: log, Additional Covariates: G83
phenotype = "G167"
transformation = log
pdfname = "G167_log_COVG83_log"

plotname = paste(phenotype, "log", sep = "_")
phenotype = transformation(as.data.frame(Pheno)[,phenotype,drop = FALSE])
addcovar = model.matrix(~sex + DOwave + diet_days + log(G83), data = Pheno)

qtl = scan1(genoprobs = genoprobs, pheno = phenotype,
            kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)

filename = paste(plot_location, pdfname, ".pdf", sep = "")
pdf(file = filename, width = 8, height = 4)
plot(qtl, map = map, main = plotname)
mtext(paste("Covariate:", "G83_log"), cex = 1, line = 0.35)
dev.off()

# Phenotype: KCl, Transformation: log, Additional Covariates: None
phenotype = "KCl"
transformation = log
pdfname = "KCl_log"

plotname = paste(phenotype, "log", sep = "_")
phenotype = transformation(as.data.frame(Pheno)[,phenotype,drop = FALSE])
addcovar = model.matrix(~sex + DOwave + diet_days, data = Pheno)

qtl = scan1(genoprobs = genoprobs, pheno = phenotype,
            kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)

filename = paste(plot_location, pdfname, ".pdf", sep = "")
pdf(file = filename, width = 8, height = 4)
plot(qtl, map = map, main = plotname)
dev.off()

# Phenotype: KCl, Transformation: log, Additional Covariates: G33
phenotype = "KCl"
transformation = log
pdfname = "KCl_log_COVG33_log"

plotname = paste(phenotype, "log", sep = "_")
phenotype = transformation(as.data.frame(Pheno)[,phenotype,drop = FALSE])
addcovar = model.matrix(~sex + DOwave + diet_days + log(G33), data = Pheno)

qtl = scan1(genoprobs = genoprobs, pheno = phenotype,
            kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)

filename = paste(plot_location, pdfname, ".pdf", sep = "")
pdf(file = filename, width = 8, height = 4)
plot(qtl, map = map, main = plotname)
mtext(paste("Covariate:", "G33_log"), cex = 1, line = 0.35)
dev.off()

# Phenotype: GLP1, Transformation: log, Additional Covariates: None
phenotype = "GLP1"
transformation = log
pdfname = "GLP1_log"

plotname = paste(phenotype, "log", sep = "_")
phenotype = transformation(as.data.frame(Pheno)[,phenotype,drop = FALSE])
addcovar = model.matrix(~sex + DOwave + diet_days, data = Pheno)

qtl = scan1(genoprobs = genoprobs, pheno = phenotype,
            kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)

filename = paste(plot_location, pdfname, ".pdf", sep = "")
pdf(file = filename, width = 8, height = 4)
plot(qtl, map = map, main = plotname)
dev.off()

# Phenotype: GLP1, Transformation: log, Additional Covariates: G83
phenotype = "GLP1"
transformation = log
pdfname = "GLP1_log_COVG83_log"

plotname = paste(phenotype, "log", sep = "_")
phenotype = transformation(as.data.frame(Pheno)[,phenotype,drop = FALSE])
addcovar = model.matrix(~sex + DOwave + diet_days + log(G83), data = Pheno)

qtl = scan1(genoprobs = genoprobs, pheno = phenotype,
            kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)

filename = paste(plot_location, pdfname, ".pdf", sep = "")
pdf(file = filename, width = 8, height = 4)
plot(qtl, map = map, main = plotname)
mtext(paste("Covariate:", "G83_log"), cex = 1, line = 0.35)
dev.off()

# Phenotype: AA, Transformation: log, Additional Covariates: None
phenotype = "AA"
transformation = log
pdfname = "AA_log"

plotname = paste(phenotype, "log", sep = "_")
phenotype = transformation(as.data.frame(Pheno)[,phenotype,drop = FALSE])
addcovar = model.matrix(~sex + DOwave + diet_days, data = Pheno)

qtl = scan1(genoprobs = genoprobs, pheno = phenotype,
            kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)

filename = paste(plot_location, pdfname, ".pdf", sep = "")
pdf(file = filename, width = 8, height = 4)
plot(qtl, map = map, main = plotname)
dev.off()

# Phenotype: AA, Transformation: log, Additional Covariates: G83
phenotype = "AA"
transformation = log
pdfname = "AA_log_COVG83_log"

plotname = paste(phenotype, "log", sep = "_")
phenotype = transformation(as.data.frame(Pheno)[,phenotype,drop = FALSE])
addcovar = model.matrix(~sex + DOwave + diet_days + log(G83), data = Pheno)

qtl = scan1(genoprobs = genoprobs, pheno = phenotype,
            kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)

filename = paste(plot_location, pdfname, ".pdf", sep = "")
pdf(file = filename, width = 8, height = 4)
plot(qtl, map = map, main = plotname)
mtext(paste("Covariate:", "G83_log"), cex = 1, line = 0.35)
dev.off()

# Phenotype: PA, Transformation: log, Additional Covariates: None
phenotype = "PA"
transformation = log
pdfname = "PA_log"

plotname = paste(phenotype, "log", sep = "_")
phenotype = transformation(as.data.frame(Pheno)[,phenotype,drop = FALSE])
addcovar = model.matrix(~sex + DOwave + diet_days, data = Pheno)

qtl = scan1(genoprobs = genoprobs, pheno = phenotype,
            kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)

filename = paste(plot_location, pdfname, ".pdf", sep = "")
pdf(file = filename, width = 8, height = 4)
plot(qtl, map = map, main = plotname)
dev.off()

# Phenotype: PA, Transformation: log, Additional Covariates: G167
phenotype = "PA"
transformation = log
pdfname = "PA_log_COVG167_log"

plotname = paste(phenotype, "log", sep = "_")
phenotype = transformation(as.data.frame(Pheno)[,phenotype,drop = FALSE])
addcovar = model.matrix(~sex + DOwave + diet_days + log(G167), data = Pheno)

qtl = scan1(genoprobs = genoprobs, pheno = phenotype,
            kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)

filename = paste(plot_location, pdfname, ".pdf", sep = "")
pdf(file = filename, width = 8, height = 4)
plot(qtl, map = map, main = plotname)
mtext(paste("Covariate:", "G167_log"), cex = 1, line = 0.35)
dev.off()

# Phenotype: PA, Transformation: log, Additional Covariates: G167 + G83
phenotype = "PA"
transformation = log
pdfname = "PA_log_COVG167_COV83_log"

plotname = paste(phenotype, "log", sep = "_")
phenotype = transformation(as.data.frame(Pheno)[,phenotype,drop = FALSE])
addcovar = model.matrix(~sex + DOwave + diet_days + log(G167) + log(G83), data = Pheno)

qtl = scan1(genoprobs = genoprobs, pheno = phenotype,
            kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)

filename = paste(plot_location, pdfname, ".pdf", sep = "")
pdf(file = filename, width = 8, height = 4)
plot(qtl, map = map, main = plotname)
mtext(paste("Covariate:", "G167_log and G83_log"), cex = 1, line = 0.35)
dev.off()

