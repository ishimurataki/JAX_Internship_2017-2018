#############################################################
# Creating summary peaks spreadsheet for ins. sec. QTL
# January 29th, 2018
# - Taki Ishimura
#############################################################
library(ggplot2)
setwd("/Users/s-ishimt/Desktop/Academic_Internship/DO_qtl")
#############################################################

# load data file
load("Data/Attie_islet_secr_data_v1.Rdata")
load("Data/pheno_clin_v6.RData")
pheno = pheno_clin
pheno_dict = pheno_clin_dict
rm(pheno_clin, pheno_clin_dict)

# setting CSV file parameters:
input.file    = "Derived_data/2018_01_24_ins_output.rds"
output.prefix = "Derived_data/2018_01_29_ins_sec"
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
qtl.summary = qtl %>% gather(key = analyte, value = lod, -chr, -pos) %>%
  group_by(analyte, chr) %>%
  filter(lod > thr) %>%
  top_n(1, lod)
qtl.summary = as.data.frame(qtl.summary[,c(3, 1:2, 4)])

write.csv(qtl.summary, paste0(output.prefix, "_qtl_summary_thresh_", thr,".csv"), 
          row.names = F, quote = F)