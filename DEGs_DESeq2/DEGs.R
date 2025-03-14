# Install the latest version of DEseq2
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")

# load the library
library(DESeq2)

#mouse data 
## Load in data
data <- read.table("datas/Mov10_full_counts.txt", header=T, row.names=1) 

meta <- read.table("meta/Mov10_full_meta.txt", header=T, row.names=1)


install.packages("pheatmap")
install.packages("DEGreport")
library(DESeq2)
library(pheatmap)
library(DEGreport)
library(tidyverse)


### Check classes of the data 
class(meta)
class(data)

#view file 
View(data)
View(meta)

#normalisation
#checking if the colnames of counts data are matched with the rownames of the meta data.. if not DESeq2 throwns an error 

### Check that sample names match in both files
all(colnames(data) %in% rownames(meta))
all(colnames(data) == rownames(meta))

#creating DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~sampletype)
# design specifies the interest of study 
View(counts(dds))


#sizefactor estimation
dds <- estimateSizeFactors(dds)

sizeFactors(dds)

#normalised counts
normalized_counts <- counts(dds, normalized=TRUE)
normalized_countss <- counts(dds)


#write normalised counts to txt file in datas directory
write.table(normalized_counts, file="datas/normalized_counts.txt", sep="\t", quote=F, col.names=NA)

## Run analysis
dds <- DESeq(dds) 
#using pre-existing size factors
#estimating dispersions
#gene-wise dispersion estimates
#mean-dispersion relationship
#final dispersion estimates
#fitting model and testing

## Check the size factors
sizeFactors(dds)


## Total number of raw counts per sample
colSums(counts(dds))


## Total number of normalized counts per sample
colSums(counts(dds, normalized=T))

#There are 3 samples namely, mov10_overexpression, mov10_knockdown, control. we compare by 1. overexpression vs control and 2. knockdown vs control 
#1. over expression vs control 
#we use contrasts to specify different sampletypes. 
library(BiocManager)
contrast_oe <- c("sampletype", "MOV10_overexpression", "control") #MOV10_overexpression compared to control condition. log2(MOV10_overespression/control)

res_tableOE_unshrunken <- results(dds, contrast = contrast_oe, alpha = 0.05)

#some genes with low mean can have higher dispersion value, in such cases, the LFC is shrunk towards 0
res_tableOE <- lfcShrink(dds, coef = "sampletype_MOV10_overexpression_vs_control", type = "apeglm")


#MA plot(mean vs log fold change)
plotMA(res_tableOE_unshrunken, ylim=c(-2,2))#unshrunken 
plotMA(res_tableOE, ylim=c(-2,2))#shrunken 


class(res_tableOE)
mcols(res_tableOE, use.names = TRUE)

res_tableOE %>% data.frame() %>% View()

#2.Knockdown vs control 
contrast_kd <- c("sampletype", "MOV10_knockdown", "control")
res_tablekd_unshrunken <- results(dds, contrast = contrast_kd, alpha = 0.05)
res_tablekd <- lfcShrink(dds, coef = "sampletype_MOV10_knockdown_vs_control", type = "apeglm")

## Summarize results
summary(res_tableOE)
summary(res_tableOE_unshrunken)

#Multiple test correction
padj.cutoff <- 0.05
lfc.cutoff <- 0.58

#table 
res_tableOE_tb <- res_tableOE %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()


str(res_tableOE_tb)# dataframe    
colnames(res_tableOE_tb) # Lists all column names

print(head(res_tableOE_tb))

head(res_tableOE_tb$padj)

#signficant genes in overexpression sample 
sigOEfiltered <- res_tableOE_tb[res_tableOE_tb$padj < padj.cutoff & 
                                  abs(res_tableOE_tb$log2FoldChange) > lfc.cutoff, ] #filter genes based on the threshold


res_tableKD_tb <- res_tablekd %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()


#signficant genes in overexpression sample 
sigKDfiltered <- res_tableKD_tb[res_tableKD_tb$padj < padj.cutoff & 
                                  abs(res_tableKD_tb$log2FoldChange) > lfc.cutoff, ]

sigOEfiltered
sigKDfiltered

# load libraries
library(ggplot2)
library(DESeq2)
library(pheatmap)

#meta and normalised counts as tibble 
library(tibble)
MOV10_meta_tb <- meta %>%
  rownames_to_column(var = "samplename") %>%
  as_tibble()

normalized_counts_tb <- normalized_counts %>% 
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

# Plot expression for single gene
plotCounts(dds, gene = "MOV10", intgroup = "sampletype")

# Save plotcounts to a data frame object
d <- plotCounts(dds, gene="MOV10", intgroup="sampletype", returnData=TRUE)

# Plotting the MOV10 normalized counts, using the samplenames 
ggplot(d, aes(x = sampletype, y = count, color = sampletype)) + 
  geom_point(position=position_jitter(w = 0.1,h = 0)) +
  geom_label(aes(label = rownames(d))) + 
  theme_bw() +
  ggtitle("MOV10") +
  theme(plot.title = element_text(hjust = 0.5))

## Plot dispersion estimates
plotDispEsts(dds)

## Order results by padj values
top20_sigOE_genes_tb <- res_tableOE_tb %>% 
  arrange(padj) %>%        # Arrange rows by padj values
  pull(gene) %>%           # Extract character vector of ordered genes
  head(n = 20) %>%         # Extract the first 20 genes
  tibble(gene = .)         # Convert to tibble and name the column "gene"

## normalized counts for top 20 significant genes
top20_sigOE_norm_tb <- normalized_counts_tb %>%
  filter(gene %in% top20_sigOE_genes_tb$gene)%>%
  as_tibble()

#gather all the sample counts for all the genes in column for plotting 
library(tidyr)
gather_top20_sigOE <- top20_sigOE_norm_tb %>%
  gather(colnames(top20_sigOE_norm_tb)[2:9], key = "samplename", value = "normalized_counts")
  
#we want our counts colored by sample group,we need to combine the metadata information with the normalized counts data into the same data frame for input to ggplot():
gathered_top20_sigOE <- inner_join(MOV10_meta_tb, gather_top20_sigOE) #the inner_join() will join both the MOV_10_meta data and top20_sig0E genes 

library(ggplot2)
## plot using ggplot2
ggplot(gathered_top20_sigOE) +
  geom_point(aes(x = gene, y = normalized_counts, color = sampletype)) +
  scale_y_log10() +
  xlab("Genes") +
  ylab("log10 Normalized Counts") +
  ggtitle("Top 20 Significant DE Genes") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5))
  
### Extract normalized expression for significant genes from the OE and control samples (4:9), and set the gene column (1) to row names
norm_OEsig <- normalized_counts_tb[,c(1,4:9)] %>% 
  filter(gene %in% sigOEfiltered$gene) %>% 
  data.frame() %>%
  column_to_rownames(var = "gene") 

### Annotate our heatmap (optional)
annotation <- MOV10_meta_tb  %>% 
  select(samplename, sampletype) %>% 
  data.frame(row.names = "samplename")

### Set a color palette
library(RColorBrewer)

heat_colors <- brewer.pal(6, "YlOrRd")

### Run pheatmap
pheatmap(norm_OEsig, 
         color = heat_colors, 
         cluster_rows = T, 
         show_rownames = F,
         annotation = annotation, 
         border_color = NA, 
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 10, 
         height = 20)

## Obtain logical vector where TRUE values denote padj values < 0.05 and fold change > 1.5 in either direction

res_tableOE_tb <- res_tableOE_tb %>% 
  mutate(threshold_OE = padj < 0.05 & abs(log2FoldChange) >= 0.58)

## Volcano plot
ggplot(res_tableOE_tb) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold_OE)) +
  ggtitle("Mov10 overexpression") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  #scale_y_continuous(limits = c(0,50)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))  



