# September 29, 2017: Program to calculate the adjusted p-values of the in vitro pancreatic insulin secretion phenotypes regressed with 
# gene expression traits. Find best correlated genes to the phenotypes of interest. 

# load all Diversity Outbred Data
DO_data <- load("/Users/s-ishimt/Desktop/Data_drive/DO_QTL/DO378_islet.RData")
DO_pheno_data <- read.csv("/Users/s-ishimt/Desktop/Data_drive/DO_QTL/usable_clin_pheno_data.csv", as.is = TRUE)

# Create dataframe of just in vitro pancreatic insulin secretion phenotypes.
Islet_ins_sec_df <- DO_pheno_data[,setdiff(intersect(grep("ins_secrete$", colnames(DO_pheno_data)), 
                                                     grep("^G", colnames(DO_pheno_data))), 
                                           grep("^GLP1_", colnames(DO_pheno_data)))]

Islet_ins_sec_df <- cbind(DO_pheno_data[,1:2],Islet_ins_sec_df)

# check which data columns need normalization
quartz()
par(mfrow=c(2,4))
for(i in 1:8){
  hist(Islet_ins_sec_df[,i+2], breaks = 20, main = colnames(Islet_ins_sec_df)[i+2], xlab =  colnames(Islet_ins_sec_df)[i+2])
}

quartz()
par(mfrow=c(2,4))
for(i in 1:8){
  hist(log(Islet_ins_sec_df[,i+2]), breaks = 20, main = colnames(Islet_ins_sec_df)[i+2], xlab =  colnames(Islet_ins_sec_df)[i+2])
} #all data columns need normalization

# log transform the clinical phenotypes to normalize the data
Islet_ins_sec_df[,3:ncol(Islet_ins_sec_df)] <- log(Islet_ins_sec_df[,3:10])

# function to compute adjusted pvalues adjusted for one covariate
pvalue_1cov <- function(x,f,c){
  fit1 <- lm(x ~ c+f)
  fit0 <- lm(x ~ c)
  anova(fit0, fit1)[2,6]
}

# create dataframe
df1 <- vector('numeric')
pvalue_1cov_allgene <- function(Islet_ins_sec_df_number){
  for(i in 1:21771){
    print(i)
    df1 <- c(df1, pvalue_1cov(Islet_ins_sec_df[,Islet_ins_sec_df_number], rankz.mrna[,i], DO_pheno_data$sex))
  }
  return(p.adjust(df1, method="BH"))
}

df2 <- data.frame(matrix(ncol = 0, nrow = 21771), stringsAsFactors = FALSE)
for(i in 1:8){
  df2 <- cbind(df2, pvalue_1cov_allgene(i+2))
  colnames(df2)[i] <- colnames(Islet_ins_sec_df)[i+2]
}
df2 <- cbind(annot.mrna$id, annot.mrna$symbol, df2)

# exportation of newly created csv file
write.csv(x=df2, file="/Users/s-ishimt/Desktop/Data_drive/DO_QTL/Pvaluetable_ins_sec.csv", row.names = FALSE)


