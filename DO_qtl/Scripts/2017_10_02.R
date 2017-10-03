# October 2-3, 2017: Program to find gene expression traits highly associated with the G33 insulin secretion phenotype, 
# then create a summative qtl for all the expression traits. A high frequency of high LOD scores at the same genomic location should
# show up on this summative qtl map. 

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
  if(pvalue_cor_test_df$G33_ins_secrete[i] > .000001){
    pvalue_cor_test_df[i, 4] <- "0"
  }
  if(pvalue_cor_test_df$G33_ins_secrete[i] <= .000001){
    pvalue_cor_test_df[i, 4] <- "1"
  }
}
colnames(pvalue_cor_test_df)[4] <- "colors"

library(ggplot2)
ggplot(data = pvalue_cor_test_df, mapping = aes(x = cor_vector, y = G33_ins_secrete)) + geom_point(aes(color = colors)) 

# identifying genes with p-values below 0.0001 and creating expression dataframe
G33_related_genes <- as.character(ins_sec_pvalues$annot.mrna.symbol[ins_sec_pvalues$G33_ins_secrete < 0.000001])
G33_related_genes_df <- cbind(DO_pheno_data[,1:2], rankz.mrna[,colnames(rankz.mrna) %in% 
                                     annot.mrna$id[annot.mrna$symbol %in% G33_related_genes]])
colnames(G33_related_genes_df)[3:ncol(G33_related_genes_df)] <- annot.mrna$symbol[annot.mrna$id %in% 
                                                                                    colnames(G33_related_genes_df)]

# running qtl analysis on all identified genes expression.
# establishing qtl2 components
library(qtl2)
library(qtl2convert)
library(ggplot2)
probs <- probs_doqtl_to_qtl2(probs = genoprobs, map = snps, chr_column = "chr", pos_column = "bp")
snps$chr <- as.character(snps$chr)
snps$chr[snps$chr=="X"] <- "20"
map <- map_df_to_list(map = snps, chr_column = "chr", pos_column = "bp")
annot.samples$diet_days <- DO_pheno_data$diet_days
K = calc_kinship(probs = probs, type = "loco", cores = 4)
addcovar <- model.matrix(~Sex + Generation + diet_days, data=annot.samples)[,-1]

# Initializing first scan as the base scan
G33_genesum_qtl <- scan1(genoprobs = probs, pheno = G33_related_genes_df[,3,drop = FALSE],
                         kinship = K, addcovar = addcovar, cores = 4, reml = TRUE)
for(z in 1:length(G33_genesum_qtl))
{
  print(z)
  if(G33_genesum_qtl[z] < 6){
    G33_genesum_qtl[z] = 0
  }
  if(G33_genesum_qtl[z] >= 6){
    G33_genesum_qtl[z] = 1
    # threshold is set to 6
  }
}

# adding scans on top of the base scan 
for(i in 4:ncol(G33_related_genes_df))
{
  print(i)
  G33_gene_qtl <- scan1(genoprobs = probs, pheno = G33_related_genes_df[,i,drop = FALSE],
                        kinship = K, addcovar = addcovar, cores = 4, reml = TRUE)
  for(z in 1:length(G33_gene_qtl))
  {
    if(G33_gene_qtl[z] < 6){
      G33_gene_qtl[z] = 0
    }
    if(G33_gene_qtl[z] >= 6){
      G33_gene_qtl[z] = 1
    }
  }
  G33_genesum_qtl <- G33_genesum_qtl + G33_gene_qtl
}

# plot the scan
plot_scan1(x = G33_genesum_qtl, map = map, lodcolumn = 1, main = "Associated gene expression traits for G33 Ins. Sec. | Threshold = 6")


