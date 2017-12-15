#call libraries
library(qtl2)
library(qtl2convert)

#Read data into working environment.
phenos37_ins_sec_data <- read.csv("/Users/s-ishimt/Desktop/Academic_Internship/DO_qtl/Data/37phenos_ins_secretion.csv", 
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

#Performing QTL analysis on all 37 clinical phenotypes
qtl <- scan1(genoprobs = probs, pheno = phenos37_ins_sec_data[,7:43,drop = FALSE],
             kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)

#Reading in qtl summary of all peaks greater than 6 LOD threshold.
phenos37_qtl_peaks_summary <- read.csv("/Users/s-ishimt/Desktop/Academic_Internship/DO_qtl/Derived_data/11_13_2017/37phenos_6LODthreshold_QTL_qtl_summary_thresh_6.csv",
                                       row.names = NULL, as.is = TRUE)

#Cutting down qtl summary to only include peaks greater than 7.5 LOD threshold.
phenos37_qtl_peaks_summary <- phenos37_qtl_peaks_summary[which(phenos37_qtl_peaks_summary$lod >7.5),]

for(i in 1:9){
  chr = phenos37_qtl_peaks_summary$chr[i]
  pheno = phenos37_qtl_peaks_summary$analyte[i]
  column = which(colnames(phenos37_ins_sec_data) == pheno)
  qtl.blup <- scan1blup(genoprobs = probs[,chr], pheno = phenos37_ins_sec_data[,column,drop = FALSE],
                        kinship = K[[chr]], addcovar = addcovar, cores = 4)
  filename = paste("/Users/s-ishimt/Desktop/Academic_Internship/DO_qtl/Plots/2017_11_27/", 
                   pheno, "_chr", chr, "_EFFECT.pdf", sep = "")
  pdf(file = filename, width = 20, height = 12)
  plot(x = qtl.blup, map = map[[chr]], columns = 1:8, col = CCcolors,
       scan1_output = qtl[,column-6, drop = FALSE], main = colnames(qtl)[column-6])
  dev.off()
}



