# load all Diversity Outbred Data
DO_data <- load("/Users/s-ishimt/Desktop/Data_drive/DO_QTL/DO378_islet.RData")
DO_pheno_data <- read.csv("/Users/s-ishimt/Desktop/Data_drive/DO_QTL/usable_clin_pheno_data.csv", as.is = TRUE)

# read in generated p-value table for insulin secretion and gene expression traits.
ins_sec_pvalues <- read.csv("/Users/s-ishimt/Desktop/Data_drive/DO_QTL/Pvaluetable_ins_sec.csv")

# find correlation of G33 insulin secretion trait with all gene expression traits, and store in vector
gene.expression.rankz <- function(gene.name){
  return(rankz.mrna[,which(colnames(rankz.mrna)==annot.mrna$id[which(annot.mrna$symbol == gene.name)[1]])])
}

cor_vector <- vector("numeric")
for(i in 1:length(annot.mrna$symbol)){
  print(i)
  gene.name <- annot.mrna$symbol[i]
  g33_ins_sec_gene_cor <- cor(gene.expression.rankz(gene.name), 
                              log(DO_pheno_data$G33_ins_secrete), use = "complete.obs")
  cor_vector <- c(cor_vector, g33_ins_sec_gene_cor)
}

# making dataframe, plotting correlation vs. p-value
pvalue_cor_test_df <- cbind(ins_sec_pvalues[,2:3], cor_vector)
for(i in 1:nrow(pvalue_cor_test_df)){
  print(i)
  if(pvalue_cor_test_df$G33_ins_secrete[i] > 0.001){
    pvalue_cor_test_df[i, 4] <- "0"
  }
  if(pvalue_cor_test_df$G33_ins_secrete[i] <= 0.001){
    pvalue_cor_test_df[i, 4] <- "1"
  }
}
colnames(pvalue_cor_test_df)[4] <- "colors"

library(ggplot2)
ggplot(data = pvalue_cor_test_df, mapping = aes(x = cor_vector, y = G33_ins_secrete)) + geom_point(aes(color = colors)) 

