rm(list = ls())

#Read data into working environment.
phenos37_ins_sec_data = read.csv("/Users/s-ishimt/Desktop/Academic_Internship/DO_qtl/Data/37phenos_MK_version2.csv", 
                                  na.strings = c("","NA"), row.names = NULL, as.is = TRUE)
rownames(phenos37_ins_sec_data) <- phenos37_ins_sec_data$mouse_id

#rankz transform function
rz.transform<-function(y) {
  rankY=rank(y, ties.method="average", na.last="keep")
  rzT=qnorm(rankY/(length(na.exclude(rankY))+1))  
  rzT
}
#rankz transform all the phenotypes
for(i in 1:37){
  phenos37_ins_sec_data[,i+6] <- rz.transform(phenos37_ins_sec_data[,i+6])
}

DO_data <- load("/Users/s-ishimt/Desktop/Data_drive/DO_QTL/DO378_islet.RData")
diet_days <- read.csv("/Users/s-ishimt/Desktop/Data_drive/DO_QTL/pheno_clin.csv", as.is = TRUE)

#retrieve diet days from other phenotype dataset
diet_days$mouse[nchar(diet_days$mouse)==5] <- sub(pattern = "DO-", replacement = "DO-0", x = diet_days$mouse[nchar(diet_days$mouse)==5])
colnames(diet_days)[1] <- "Mouse.ID"
diet_days <- diet_days[match(x=annot.samples$Mouse.ID, table=diet_days$Mouse.ID),]
stopifnot(annot.samples$Mouse.ID == diet_days$Mouse.ID)
diet_days <- diet_days[,which(colnames(diet_days) %in% c("Mouse.ID", "diet_days"))]

#calculate qtl mapping elements
probs <- probs_doqtl_to_qtl2(probs = genoprobs, map = snps, chr_column = "chr", pos_column = "bp")
snps$chr <- as.character(snps$chr)
snps$chr[snps$chr=="X"] <- "20"
map <- map_df_to_list(map = snps, chr_column = "chr", pos_column = "bp")
annot.samples$diet_days <- diet_days$diet_days
K = calc_kinship(probs = probs, type = "loco", cores = 4)

addcovar <- model.matrix(~Sex + Generation + diet_days, data=annot.samples)[,-1]

#Reading in qtl summary of all peaks greater than 6 LOD threshold.
phenos37_qtl_peaks_summary <- read.csv("/Users/s-ishimt/Desktop/Academic_Internship/DO_qtl/Derived_data/11_13_2017/37phenos_6LODthreshold_summary.csv",
                                       row.names = NULL, as.is = TRUE)

#Cutting down qtl summary to only include peaks greater than 7.5 LOD threshold.
phenos37_qtl_peaks_summary <- phenos37_qtl_peaks_summary[which(phenos37_qtl_peaks_summary$lod >7.5),]

#source mediation analysis functions
source("/Users/s-ishimt/Desktop/Academic_Internship/DO_qtl/Scripts/Mediation_Analysis_Intermediate.R")

#set up different function arguments.
mediator = rankz.mrna
annotation = annot.mrna[,-c(4,5,6,8,9)]
colnames(annotation)[4] = "pos"
annotation[,4] = annotation[,4]/1e+06
for(i in 1:nrow(annotation)){
  if(annotation$chr[i] == "MT"){
    annotation$chr[i] = "M"
  }
} #Renames all mitochondrial genes labeled as MT as M for function readability. 
covar = addcovar
load("/Users/s-ishimt/Desktop/Data_drive/DO_QTL/mouse.chrlen.rda")

for(i in 1:nrow(phenos37_qtl_peaks_summary)){
  phenotype = phenos37_qtl_peaks_summary$analyte[i]
  target = as.matrix(phenos37_ins_sec_data[,colnames(phenos37_ins_sec_data) == phenotype, drop = FALSE])
  geno_label = paste(as.character(phenos37_qtl_peaks_summary$chr[i]), "_",
                     as.character(phenos37_qtl_peaks_summary$pos[i]*1e+06), sep = "")
  qtl.geno = genoprobs[,,dimnames(genoprobs)[[3]] == geno_label]
  med <- mediation.scan(target = target,
                        mediator = mediator,
                        annotation = annotation,
                        covar = covar,
                        qtl.geno = qtl.geno)
  quartz()
  plot.mediation(med)
}


