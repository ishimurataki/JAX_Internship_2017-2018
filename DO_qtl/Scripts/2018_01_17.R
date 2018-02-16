#############################################################
# Insulin secretion data general analysis
# January 17th, 2018
# - Taki Ishimura
#############################################################
library(ggplot2)
setwd("/Users/s-ishimt/Desktop/Academic_Internship/DO_qtl")
#############################################################

# load data file
load("Data/Attie_islet_secr_data_v1.Rdata")
<<<<<<< HEAD
load("Data/pheno_clin_v6.RData")
pheno = pheno_clin
pheno_dict = pheno_clin_dict
rm(pheno_clin, pheno_clin_dict)
=======
>>>>>>> cf59327ef1884d2eae4729d9e8d5d2772e79e36c

# initialize phenotypes vector
phenotypes = c("num_islets",
               "Ins_per_islet",
               "WPIC",
               "Glu0",
               "Ins0",
               "HOMA_B",
               "HOMA_IR",
               "Glu_tAUC",
               "Glu_iAUC",
               "Ins_tAUC",
               "Ins_iAUC",
               "Glu",
               "Ins",
               "TG",
               "Gcg_content",
               "Gcg_secreted",
               "G33",
               "G83",
               "G167",
               "KCl",
               "GLP1",
               "AA",
               "PA")

# re-name the phenotypes
renaming <- read.csv("QTL_parameters/18_01_02_ins_renaming_parameters.csv", stringsAsFactors=FALSE)
phe_col <- colnames(pheno)
phe_col[match(renaming$old_name, phe_col)] <- renaming$new_name
colnames(pheno) <- phe_col

# Q-Q plots of the phenotypes and save in pdf
pdf(file = "Plots/2018_01_17/ins_sec_qqplots.pdf", width = 5, height = 5)
plot = qplot(sample = sqrt(pheno[,phenotypes[1]]), main = phenotypes[1], color = "cy1")
plot(plot)
for(i in 2:length(phenotypes))
{
  pheno_data = log(pheno[,phenotypes[i]])
  plot = qplot(sample = pheno_data, main = phenotypes[i],  color = "cy1")
  plot(plot)
}
dev.off()

correlation =  cor(log(pheno$Ins0), log(pheno$Glu0), use = "complete.obs")
plot = ggplot(data = pheno, aes(x = log(Ins0), y = log(Glu0), color = sex, group = sex)) +
  geom_smooth(method = 'lm', formula = y~x) +
  geom_point() +
  labs(x = "Ins0", y = "Glu0", title = "Glu0 vs Ins0", subtitle = correlation)

# read parameters file for creating scatterplots
scatterplots = read.csv("Plots/2018_01_17/scatterplots.csv", stringsAsFactors = FALSE)
pdf(file = "Plots/2018_01_17/ins_sec_scatterplots.pdf", width = 5, height = 5)
for(i in 1:nrow(scatterplots))
{
<<<<<<< HEAD
  quartz()
=======
>>>>>>> cf59327ef1884d2eae4729d9e8d5d2772e79e36c
  xname = scatterplots$X[i]
  yname = scatterplots$Y[i]
  x = log(pheno[,xname])
  y = log(pheno[,yname])
  data = data.frame(x, y, sex = pheno$sex)
  correlation = cor(x,y, use = "complete.obs")

  plot = ggplot(data = data,aes(x = x, y = y, color = sex, group = sex)) +
    geom_smooth(method = 'lm', formula = y~x) + 
    geom_point() +
    labs(x = xname, y = yname, title = paste(yname, "vs", xname, sep = " "), subtitle = correlation)
  
  plot(plot)
}
dev.off()

