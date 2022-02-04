# figure 2 expression
library(stringr)
library(reshape2)
library(tidyverse)
expression <- read.csv("hv_R3_genes_norm_counts_de_tax_sp.csv", stringsAsFactors = F)
expression <- expression[grep("HORVU", expression$gene, ignore.case=TRUE), ]

gi = expression[,2:91]
gi = log(expression[,2:91]+1)

rownames(gi)<- expression$gene

# Calculating means for gi
gmeans = matrix(nrow = nrow(gi), ncol=30)
gmeans<- data.frame(gmeans)
j = 1
for(k in 1:30){
  gmeans[,k] = rowMeans((gi[,j:(j+2)]))
  colnames(gmeans)[k]<- str_sub(colnames(gi)[j], end=-4)
  j = j+3
}
rm(j, k, gi)

# Reading in means for genotype, depending on the order in the table

# For V2 new annotation 
means.bln1.norm = gmeans[,19:24]
means.dm.norm = gmeans[,25:30]
means.mla6.norm = gmeans[,7:12]
means.rar3.norm = gmeans[,13:18]
means.wt.norm = gmeans[,1:6]

gmeans$gene<- expression$gene
gmeans_long<- melt(gmeans,id.vars = "gene")
gmeans_long$genotype<- sapply(gmeans_long$variable, FUN=function(x)str_split(x,"_")[[1]][1])
gmeans_long$timepoint<- sapply(gmeans_long$variable, FUN=function(x)str_split(x,"_")[[1]][2])
gmeans_long$WT_value<- rep(gmeans_long$value[gmeans_long$genotype=="wt"], 5)

#str_sub(gmeans_long$variable, end=-4)
genotype_labels<- c("bln1", "mla6/bln1", "mla6", "sgt1")
names(genotype_labels)<- c("bln1", "dm", "mla6", "rar3")
timepoint_labels<-c("0 HAI", "16 HAI", "20 HAI", "24 HAI", "32 HAI", "48 HAI")
names(timepoint_labels)<- c("t0", "t16", "t20", "t24", "t32", "t48")


#calculate correlations

gmean_long_mutants<- gmeans_long[gmeans_long$genotype!="wt",]
correlations<- gmean_long_mutants %>% group_by(genotype, timepoint) %>% summarise(cor=cor(value, WT_value))
correlations$cor<- round(correlations$cor, digits = 4)

library(tidyverse)
pdf("~/iowa_state/lab/MLA GRN publication/GRN paper/GRN_analysis/pairwise_expression.pdf", width = 18, height = 12, fonts = "ArialMT", pointsize = 18)
# ggplot(gmeans_long[gmeans_long$genotype!="wt",], aes(x=WT_value, y=value)) +geom_point(size=0.5) + 
#   facet_grid(genotype~timepoint, labeller = labeller(genotype = genotype_labels, timepoint = timepoint_labels)) + 
#   theme_bw(base_size = 40) + theme(strip.text.y = element_text(face = "bold.italic"), strip.text.x = element_text(face = "bold")) + 
#   xlab("Log(WT-expression)") + ylab("Log(Mutant-expression)")

ggplot(gmeans_long[gmeans_long$genotype!="wt",], aes(x=WT_value, y=value)) +geom_point(size=0.5) + 
  facet_grid(genotype~timepoint, labeller = labeller(genotype = genotype_labels, timepoint = timepoint_labels)) + 
  theme_bw(base_size = 40) + theme(strip.text.y = element_text(face = "bold.italic"), strip.text.x = element_text(face = "bold")) + 
  xlab("Log(WT-expression)") + ylab("Log(Mutant-expression)") + geom_text(x = 3.5, y = 12,data=correlations, aes(label = cor), size=9)

dev.off()


