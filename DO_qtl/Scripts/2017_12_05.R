###################################################################################
# 12/05/2017 - 12/10/2017
# Program to quickly run through many QTL scans with varying parameters.
# Handles different phenotypes to be mapped, differed data normalization methods,
# different covariates to be used in scans, and scan grouping in the exported pdf.
# - Taki Ishimura
###################################################################################

rm(list = ls())

library(qtl2)
library(qtl2convert)

#alright, time to get this party started.
#load in all the phenotype data from the Attie Lab.
phenodata = load("/Users/s-ishimt/Desktop/Academic_Internship/DO_qtl/Data/pheno_clin_v5.RData")
probs = readRDS("/Users/s-ishimt/Desktop/Academic_Internship/DO_qtl/Data/attie_DO500_genoprobs_v5.rds")
map = readRDS("/Users/s-ishimt/Desktop/Academic_Internship/DO_qtl/Data/grid_pmap.rds")

#Read in QTL mapping parameters
QTL_parameters = read.csv("/Users/s-ishimt/Desktop/Academic_Internship/DO_qtl/QTL_parameters/17_12_05_parameters.csv", 
                          stringsAsFactors = FALSE)

#Tailoring QTL phenotype parameters to fit the input data.
QTL_parameters[QTL_parameters == "Glu0"]             = "Glu_0min"
QTL_parameters[QTL_parameters == "Ins0"]             = "Ins_0min"
QTL_parameters[QTL_parameters == "HOMA_B"]           = "HOMA_B_0min"
QTL_parameters[QTL_parameters == "HOMA_IR"]          = "HOMA_IR_0min"
QTL_parameters[QTL_parameters == "Glu_tAUC"]         = "Glu_tAUC"
QTL_parameters[QTL_parameters == "Glu_iAUC"]         = "Glu_iAUC"
QTL_parameters[QTL_parameters == "Ins_tAUC"]         = "Ins_tAUC"
QTL_parameters[QTL_parameters == "Ins_iAUC"]         = "Ins_iAUC"
QTL_parameters[QTL_parameters == "G33_ins"]          = "G33_ins_secrete"
QTL_parameters[QTL_parameters == "G33_fract"]        = "G33_fract_ins_secrete"
QTL_parameters[QTL_parameters == "G83_ins"]          = "G83_ins_secrete"
QTL_parameters[QTL_parameters == "G83_FC"]           = "G83_FC_ins_secrete"
QTL_parameters[QTL_parameters == "G83_fract"]        = "G83_fract_ins_secrete"
QTL_parameters[QTL_parameters == "G167_ins"]         = "G167_ins_secrete"
QTL_parameters[QTL_parameters == "G167_FC"]          = "G167_FC_ins_secrete"
QTL_parameters[QTL_parameters == "G167_fract"]       = "G167_fract_ins_secrete"
QTL_parameters[QTL_parameters == "KCl_ins"]          = "KCl_G33_ins_secrete"
QTL_parameters[QTL_parameters == "KCl_FC"]           = "KCl_G33_FC_ins_secrete"
QTL_parameters[QTL_parameters == "KCl_fract"]        = "KCl_G33_fract_ins_secrete"
QTL_parameters[QTL_parameters == "GLP1_ins"]         = "GLP1_G83_ins_secrete"
QTL_parameters[QTL_parameters == "GLP1_FC_G83"]      = "GLP1_G83_FC_G83_ins_secrete"
QTL_parameters[QTL_parameters == "GLP1_fract"]       = "GLP1_G83_fract_ins_secrete"
QTL_parameters[QTL_parameters == "AA_ins"]           = "AA_G83_ins_secrete"
QTL_parameters[QTL_parameters == "AA_FC_G83"]        = "AA_G83_FC_G83_ins_secrete"
QTL_parameters[QTL_parameters == "AA_fract"]         = "AA_G83_fract_ins_secrete"
QTL_parameters[QTL_parameters == "PA_ins"]           = "PA_G167_ins_secrete"
QTL_parameters[QTL_parameters == "PA_FC_G167"]       = "PA_G167_FC_G167_ins_secrete"
QTL_parameters[QTL_parameters == "PA_fract"]         = "PA_G167_fract_ins_secrete"
QTL_parameters[QTL_parameters == "Gcg_content"]      = "Gcg_content"
QTL_parameters[QTL_parameters == "Gcg"]              = "Gcg_secreted"
QTL_parameters[QTL_parameters == "Gcg_percent"]      = "Gcg_secreted_percent"
QTL_parameters[QTL_parameters == "num_islets"]       = "num_islets"
QTL_parameters[QTL_parameters == "Ins_per_islet"]    = "Ins_per_islet"
QTL_parameters[QTL_parameters == "WPIC"]             = "WPIC"
QTL_parameters[QTL_parameters == "Glu"]              = "Glu_sac"
QTL_parameters[QTL_parameters == "Ins"]              = "Ins_sac"
QTL_parameters[QTL_parameters == "TG"]               = "TG_sac"

#checking the parameters to see whether they are acceptable.
for(i in 1:nrow(QTL_parameters)){
  print(i)
  stopifnot(QTL_parameters$Phenotype[i] %in% colnames(pheno_clin))
  stopifnot(QTL_parameters$trans[i] %in% c("log", "rankZ"))
  stopifnot(QTL_parameters$covariate[i] %in% c(colnames(pheno_clin),NA))
  stopifnot(QTL_parameters$cov_trans[i] %in% c("log", "rankZ", NA))
}

#rankz transform function
rz.transform<-function(y) {
  rankY=rank(y, ties.method="average", na.last="keep")
  rzT=qnorm(rankY/(length(na.exclude(rankY))+1))  
  rzT
}

#function to replace rankz transformed data in single column phenotype data with row identification.
rz.transform.phenotype = function(y){
  y[,1] = as.numeric(rz.transform(y))
  y
}

#make sure that the order of the mice in the probability matrix is the same as the the order of the mice in the phenotype data.
for(i in 1:20){
  prob_order = which(rownames(probs[[i]]) %in% pheno_clin$mouse)
  probs[[i]] = probs[[i]][prob_order,,]
}

#calculate the kinship matrix
K = calc_kinship(probs = probs, type = "loco", cores = 4)

#create covariate matrix elements
pheno_clin$sex = as.factor(pheno_clin$sex)
pheno_clin$DOwave = as.factor(pheno_clin$DOwave)

#HERE WE GOO!! This is where it gets exciting. 
#Using entries in the QTL_parameters spreadsheet to loop through various QTL scans.
group = 1
location = "/Users/s-ishimt/Desktop/Academic_Internship/DO_qtl/Plots/2017_12_12/" 
filename = paste(location, "group", group, ".pdf", sep = "")
pdf(file = filename, width = 8, height = 4)
for(i in 1:60){
  print(i)
  
  phenotype = QTL_parameters$Phenotype[i]
  pheno_transformation = QTL_parameters$trans[i]
  covariate = QTL_parameters$covariate[i]
  cov_transformation = QTL_parameters$cov_trans[i]
  
  #normalizing the phenotype data by log or rankZ depending on the input parameter. 
  if(pheno_transformation == "log"){
    phenotype = log(pheno_clin[,phenotype,drop = FALSE])
  } else if(pheno_transformation == "rankZ"){
      phenotype = rz.transform.phenotype(pheno_clin[,phenotype,drop = FALSE])
    }
  
  #normalizing the covariate data by log or rankZ depending on the input parameter.
  #creates matrix of covariates, always including sex, wave, and diet days
  if(is.na(cov_transformation) == FALSE){
    if(cov_transformation == "log"){
      covname = paste(covariate, "log", sep = "")
      covariate = log(pheno_clin[,covariate, drop = FALSE])
    } else if(cov_transformation == "rankZ"){
        covname = paste(covariate, "rankZ", sep = "")
        covariate = rz.transform.phenotype(pheno_clin[,covariate, drop = FALSE])
      }
    addcovar = model.matrix(~sex + DOwave + diet_days + covariate[,1], data = pheno_clin)[,-1]
  } else{
      covname = NA
      addcovar = model.matrix(~sex + DOwave + diet_days, data = pheno_clin)[,-1]
  }
  
  #running the genome scan.
  qtl = scan1(genoprobs = probs, pheno = phenotype,
               kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)
  
  #plotting qtl in pdf.
  plot(qtl, map = map, main = paste(colnames(qtl), pheno_transformation, sep = ""))
  if(is.na(covname) == FALSE){
    mtext(paste("Covariate:", covname), cex = 1, line = 0.35)
  } #inserting covariate name in plot if a covariate exists.
 
  #clustering the qtl plots on different pdfs.
  #depends on the "breaks" parameter in the QTL_parameters spreadsheet. 
  Newplot = QTL_parameters$break.[i]
  if(Newplot){
    dev.off()
    group = group + 1
    filename = paste(location, "group", group, ".pdf", sep = "")
    pdf(file = filename, width = 8, height = 4)
  }
  
  if(is.na(covname) == FALSE){
    colnames(qtl) = paste(colnames(qtl), pheno_transformation, "_COV", covname, sep = "")
  } else{
    colnames(qtl) = paste(colnames(qtl), pheno_transformation, sep = "")
  }
  
  #initializing datatable of QTL values if the loop is on its first cycle.
  #Otherwise, appending QTL values of succesive scans to the first datafile.
  if(i == 1){
    qtldatatable = qtl
  } else{
    qtldatatable = cbind(qtldatatable, qtl)
  }
}
dev.off()
saveRDS(object = qtldatatable, 
          file = "/Users/s-ishimt/Desktop/Academic_Internship/DO_qtl/Derived_data/2017_12_12/qtldatatable.rds")
