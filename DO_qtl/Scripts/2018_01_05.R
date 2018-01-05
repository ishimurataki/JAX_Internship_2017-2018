#############################################################
# High LOD peaks for insulin secretion phenotypes QTL
# January 5th, 2018
# - Taki Ishimura
#############################################################
library(qtl2)
library(qtl2convert)
library(tidyverse)
setwd("/Users/s-ishimt/Desktop/Academic_Internship/DO_qtl")
#############################################################
# load data file
load("Data/Attie_islet_secr_data_v1.Rdata")
ls()

# setting CSV file parameters:
input.file    = "Derived_Data/2018_01_02_ins_output.rds"
output.prefix = "Derived_data/2018_01_05_ins_sec"
thr           = 7

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

## Plotting chr 11
plot_list = c("Ins0_log", "Ins0_log_covlogGlu0m50", "Ins_tAUC_log", "HOMA_B_log", "Ins0_log_covlogGlu0")
quartz()
par(mfrow = c(3,2))
for(i in 1:length(plot_list)){
  plot(qtl, lodcolumn = plot_list[i], map = map, main = plot_list[i], chr = 11)
}

cor(log(pheno$Ins_0min), log(pheno$Ins_tAUC), use = "complete.obs")



