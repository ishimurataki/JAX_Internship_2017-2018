################################################################################
# 37 Phenotypes from Mark Keller: Data Importation, Clean Up, and normalization#
################################################################################

#Read data into working environment. Change all empty or NA cells into NAs.
phenos37_ins_sec_data <- read.csv("/Users/s-ishimt/Desktop/Academic_Internship/DO_qtl/Data/37phenos_ins_secretion.csv", 
                        na.strings = c("","NA"), row.names = NULL)
View(phenos37_ins_sec_data) #view data
phenos37_ins_sec_data$mouse_id <- as.character(phenos37_ins_sec_data$mouse_id)

#Read other DO data in working environment.
DO_data <- load("/Users/s-ishimt/Desktop/Data_drive/DO_QTL/DO378_islet.RData")

#subset new imported phenotype data so that individuals match the expression data
phenos37_ins_sec_data <- 
  phenos37_ins_sec_data[c(which(phenos37_ins_sec_data$mouse_id %in% annot.samples$Mouse.I)),]

# sort 37 phenotype data so that the order of individuals line up with 
phenos37_ins_sec_data <- 
  phenos37_ins_sec_data[c(order(phenos37_ins_sec_data$animal_.)),]
identical(annot.samples$Mouse.ID, phenos37_ins_sec_data$mouse_id)

#save edited file to desktop directory
write.csv(phenos37_ins_sec_data, file = "/Users/s-ishimt/Desktop/Academic_Internship/DO_qtl/Data/37phenos_ins_secretion.csv", row.names = FALSE)

#check normality of 37 phenotypes data
quartz()
par(mfrow = c(5,8))
for(i in 1:37)
{
  pheno <- log(phenos37_ins_sec_data[,i+6])
  hist(pheno)
}

#loop to check which phenotype values are negative
for(i in 1:37){
  pheno <- phenos37_ins_sec_data[,i + 6]
  for(z in 1:378){
    if(is.na(pheno[z]) == FALSE & pheno[z] < 0){
      print("*************")
      print(colnames(phenos37_ins_sec_data)[i+6]) 
      print(z)
      print(phenos37_ins_sec_data[z,i+6])
    }
  }
} #single sample in HOMA_B and multiple samples in Ins_iAUC are negative

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

#check normalization:
quartz()
par(mfrow = c(5,8))
for(i in 1:37)
{
  pheno <- phenos37_ins_sec_data[,i+6]
  hist(pheno)
}
