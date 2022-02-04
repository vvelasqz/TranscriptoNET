# Code by Valeria velasquez to analyse the data using the barley V3 and bgh v2 annotation and taxon specific normalization 
#https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html

# Load libraries
library(DESeq2)
library(gtools)

# Read, select, group and format count data
countData=read.table("HvR3_Bgh_merged_salmon_gene_counts.tsv", header=TRUE, row.names = 1)[,-1]
colnames(countData)
# Round count matrix
countData[]=round(as.matrix(countData))

# Selecting genes by counts
counts = subset(countData, rowSums(countData > 1) >= 3)
print(nrow(counts))
#rm(countData)

# Creating sample info and reordering 
time=c(0, 16, 20, 24, 32, 48)
gen=c("wt", "mla6", "rar3","bln1", "dm")
counts<- counts[, c(grep("wt", colnames(counts)),
                    grep("mla6", colnames(counts)),
                    grep("rar3", colnames(counts)),
                    grep("bln1", colnames(counts)),
                    grep("dm", colnames(counts)))]
colnames(counts)
cols=data.frame(genotype=rep(gen, each=3*6), timepoint=rep(rep(time, each=3),5), replication=rep(c(1,2,3),5*6))
cols$timepoint = as.factor(cols$timepoint)
cols$replication = as.factor(cols$replication)
cols$group <- factor(paste0(cols$genotype,"_t", cols$timepoint))
colnames(counts) <- paste0(rep(gen, each=3*6), "_t",rep(rep(time, each=3),5), "_r",rep(c(1,2,3),5*6))
rownames(cols) = colnames(counts)

############################################################################################
# Creating DESeqDataSet and running DESeq for_de counts for each taxon https://support.bioconductor.org/p/66067/
taxon = list(c("BLGH"), 
             c("HORVU.MOREX.r3.1", "HORVU.MOREX.r3.2", "HORVU.MOREX.r3.3", "HORVU.MOREX.r3.4", "HORVU.MOREX.r3.5", 
               "HORVU.MOREX.r3.6", "HORVU.MOREX.r3.7", "HORVU.MOREX.r3.Un"))

counts_taxsp<-counts
counts_taxsp[!is.na(counts_taxsp)]<- NA
#counts_taxsp_na<- counts_taxsp
for (sp in 1:2) {
  subtable<- subset(counts, rownames(counts) %in% grep(paste(taxon[[sp]], collapse="|"), rownames(counts), value = TRUE)) 
  dds_de=DESeqDataSetFromMatrix(countData = subtable, colData = cols, design = ~ group)
  dds_de=estimateSizeFactors(dds_de)
  #vsd_de <- vst(dds_de)
  counts_taxsp[rownames(counts) %in% grep(paste(taxon[[sp]], collapse="|"), rownames(counts), value = TRUE),] <- subtable/rep(sizeFactors(dds_de), each = (nrow(subtable)))
  #counts_taxsp_na[rownames(counts) %in% grep(paste(taxon[[sp]], collapse="|"), rownames(counts), value = TRUE),] <- assay(vsd_de)
}

counts_taxsp<- round(counts_taxsp)
#counts_taxsp_na <- round(counts_taxsp_na)
# Extracting and writing normalized counts
norm_counts_de = data.frame(gene = rownames(counts_taxsp), counts_taxsp, row.names = NULL)
norm_counts_de[norm_counts_de$gene=="BLGH_06491",]
#norm_counts_na = data.frame(gene = rownames(counts_taxsp_na), counts_taxsp_na, row.names = NULL)

write.csv(norm_counts_de, "hv_R3_genes_norm_counts_de_tax_sp.csv", row.names = FALSE)

# now variance stabilized data for network analysis or clustering
#write.csv(norm_counts_na, "hv_V2_genes_norm_counts_na_tax_sp.csv", row.names = FALSE)

##############################################################################################
# Creating DESeqDataSet and running DESeq for _de per timepoint counts removed the interactions for replicates
# The design matrix has the same number of samples and coefficients to fit,
# so estimation of dispersion is not possible. Treating samples
# as replicates was deprecated in v1.20 and no longer supported since v1.22.

dds=DESeqDataSetFromMatrix(countData = counts_taxsp, colData = cols, design = ~ group)
normalizationFactors(dds) <- matrix(1,ncol = ncol(counts_taxsp), nrow = nrow(counts_taxsp))
#vsd <- round(assay(vst(dds)))
dds <- DESeq(dds, quiet = TRUE)

#Now the contrasts between genotypes and timepoints https://support.bioconductor.org/p/67600/#67612
resultsNames(dds)

dds<- readRDS("~/iowa_state/lab/RNAseq/DESeq2 HvR3/dds_tax_sp_R3.RDS")

r1<-results(dds, contrast=c("group", "wt_t16", "bln1_t16"))
r2<- r1$baseMean
# r2<-results(dds, contrast=c("group", "t20_wt", "t20_mla6"))

combi <-combn(c("wt", "mla6", "rar3","bln1", "dm"),2)
padj_table <- data.frame(gene=rownames(counts_taxsp))
fc_table <- data.frame(gene=rownames(counts_taxsp))
FC_padj_table <- data.frame(gene=rownames(counts_taxsp))
counter<-0
contrast_tables<- list()
for(t in paste0("t",c(0, 16, 20, 24, 32, 48))){
  for(i in 1:10){
    counter=counter+1
    print(counter)
    name<- paste0(t,"_", combi[1,i],"_",  combi[2,i])
    contrast_tables[[name]]<-results(dds, contrast=c("group", paste0(combi[1,i], "_",t ), paste0(combi[2,i], "_", t )))
    padj_table<- cbind(padj_table, contrast_tables[[name]]$padj)
    colnames(padj_table)[1+counter]<- paste0(name, "_padj")
    fc_table<- cbind(fc_table, contrast_tables[[name]]$log2FoldChange)
    colnames(fc_table)[1+counter]<- paste0(name, "_log2FoldChange")
    FC_padj_table<- cbind(FC_padj_table, contrast_tables[[name]][,c("log2FoldChange","padj")])
    colnames(FC_padj_table)[(2*counter):(2*counter+1)]<- c(paste0(name, "_log2FoldChange"), paste0(name, "_padj"))
  }
}

write.csv(padj_table, "hv_R3_genes_deseq2_results_pairwise_genotype_padj_tax_sp.csv", row.names = FALSE)
write.csv(fc_table, "hv_R3_genes_deseq2_results_pairwise_genotype_logfc_tax_sp.csv", row.names = FALSE)
write.csv(FC_padj_table, "hv_R3_genes_deseq2_results_pairwise_genotype_logfc_padj_tax_sp.csv", row.names = FALSE)
saveRDS(contrast_tables, "contrast_tables_pairwise_tax_sp_hv_R3.RDS")
saveRDS(dds, "dds_tax_sp_R3.RDS")


# deciding on a threshold for significant comparisons
library(tidyverse)
contrast_tables_pairwise_tax_sp<- readRDS("contrast_tables_pairwise_tax_sp_hv_R3.RDS")
chrs = c("HORVU.MOREX.r3.1", "HORVU.MOREX.r3.2", "HORVU.MOREX.r3.3", "HORVU.MOREX.r3.4", "HORVU.MOREX.r3.5", 
         "HORVU.MOREX.r3.6", "HORVU.MOREX.r3.7", "HORVU.MOREX.r3.Un")
chrs2 = c("BLGH")
thr_vals<- c(0.05, 0.01, 0.005, 0.003,0.0025, 0.001, 0.0008)
sig_tests<- data.frame(org=c(rep("Hv", 60), rep("Bgh", 60)), test=rep(names(contrast_tables_pairwise_tax_sp),2))
for (names in names(contrast_tables_pairwise_tax_sp)) {
  t0<-contrast_tables_pairwise_tax_sp[[names]]
  t0_1<- subset(t0, rownames(t0) %in% grep(paste(chrs, collapse="|"), rownames(t0), value = TRUE))
  t0_2<- subset(t0, rownames(t0) %in% grep(paste(chrs2, collapse="|"), rownames(t0), value = TRUE))
  for(thr in thr_vals){
    sig_tests[sig_tests$org=="Hv" & sig_tests$test==names, paste0("n_",thr)]<- sum(t0_1$padj<thr, na.rm = T)*thr
    sig_tests[sig_tests$org=="Bgh" & sig_tests$test==names, paste0("n_",thr)]<- sum(t0_2$padj<thr, na.rm = T)*thr
    # print(paste0("num of sig tests for barley ",names,", ",thr, ", ",sum(t0_1$padj<thr, na.rm = T)*thr))
    # print(paste0("num of sig tests for blumeria ",names,", ",thr, ", ",sum(t0_2$padj<thr, na.rm = T)*thr))
  }
}

sig_tests %>% group_by(org) %>% summarise(mean_0.05=mean(n_0.05), mean_0.01=mean(n_0.01), mean_0.005=mean(n_0.005),
                                          mean_0.003=mean(n_0.003), mean_0.0025=mean(n_0.0025), mean_0.001=mean(n_0.001),
                                          mean_0.0008=mean(`n_8e-04`))


# A tibble: 2 x 8
# org   mean_0.05 mean_0.01 mean_0.005 mean_0.003 mean_0.0025 mean_0.001 mean_0.0008
# <chr>     <dbl>     <dbl>      <dbl>      <dbl>       <dbl>      <dbl>       <dbl>
#   1 Bgh        32.3      4.19       1.80      0.973       0.782      0.266       0.205
# 2 Hv        143.      19.0        8.22      4.47        3.60       1.22        0.939

#Bgh 0.003
#hv 0.0008