
#call libraries
library(qtl2)
library(qtl2convert)
library(tidyverse)

#Read data into working environment. Change all empty or NA cells into NAs.
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

View(phenos37_ins_sec_data) #view data

#Read other DO data in working environment.
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

#Performing QTL analysis on all 37 clinical phenotypes and exporting as pdf:
for(i in 1:37){
  filename = paste("/Users/s-ishimt/Desktop/Academic_Internship/DO_qtl/Plots/2017_11_13/", 
                   colnames(phenos37_ins_sec_data)[i+6], "_LODscan.pdf", sep = "")
  qtl <- scan1(genoprobs = probs, pheno = phenos37_ins_sec_data[,i+6,drop = FALSE],
               kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)
  pdf(file = filename, width = 15, height = 6)
  plot_scan1(x = qtl, map = map, lodcolumn = 1, main = colnames(phenos37_ins_sec_data)[i+6])
  dev.off()
}

#Creating csv file of all 37 phenotype QTL scans
qtl <- scan1(genoprobs = probs, pheno = phenos37_ins_sec_data[,7:43,drop = FALSE],
             kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)
saveRDS(qtl, "/Users/s-ishimt/Desktop/Academic_Internship/DO_qtl/Derived_data/37phenos_QTL.rds")

################################################
#Creating csv file of LOD scores higher than 6 #
################################################
input.file    = "/Users/s-ishimt/Desktop/Academic_Internship/DO_qtl/Derived_data/37phenos_QTL.rds"
output.prefix = "/Users/s-ishimt/Desktop/Academic_Internship/DO_qtl/Derived_data/11_13_2017/37phenos"
thr           = 6

print(paste("INPUT FILE =", input.file))
print(paste("OUTPUT PREFIX =", output.prefix))
print(paste("THR =", thr))

# Read in the aggregated QTL file.
qtl = readRDS(file = input.file)
qtl = data.frame(qtl)

chr_pos <- do.call(rbind, (strsplit(rownames(qtl), "_")))

qtl= data.frame(chr = chr_pos[,1],
                pos = as.numeric(chr_pos[,2])/1e+06, qtl)
qtl.summary = qtl %>% gather(analyte, lod, -chr, -pos) %>%
  group_by(analyte, chr) %>%
  filter(lod > thr) %>%
  top_n(1, lod)
qtl.summary = as.data.frame(qtl.summary[,c(3, 1:2, 4)])

write.csv(qtl.summary, paste0(output.prefix, "_qtl_summary_thresh_", thr,".csv"), 
          row.names = F, quote = F)


