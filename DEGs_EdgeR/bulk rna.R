library(magrittr)

library(statmod)

library(edgeR)

library(org.Mm.eg.db)

library(tidyverse)

library(ggplot2)

targets <- read.delim("targets.txt", stringsAsFactors=FALSE)
#merge cell types and status
group <- paste(targets$CellType, targets$Status, sep=".")
group <- factor(group)
table(group)

#2 types of cells - lactating and basal stem cells, 3 types of conditions namely PREGNANT, LACTATING and VIRGIN, 2 biological replicates.. 2 * 3 * 2 = 12 samples 
FileURL <- paste("http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE60450","format=file","file=GSE60450_Lactation-GenewiseCounts.txt.gz", sep="&")
download.file(FileURL, "GSE60450_Lactation-GenewiseCounts.txt.gz")

GenewiseCounts <- read.delim("GSE60450_Lactation-GenewiseCounts.txt.gz",row.names="EntrezGeneID")

colnames(GenewiseCounts) <- substring(colnames(GenewiseCounts),1,7)# subset the colname
dim(GenewiseCounts)
head(GenewiseCounts)

library(edgeR)
y <- DGEList(GenewiseCounts[,-1], group=group, genes=GenewiseCounts[,1,drop=FALSE])
options(digits=3)
y$samples             


install.packages("org.Mm.eg.db")
y$genes$Symbol <- mapIds(org.Mm.eg.db, rownames(y),keytype="ENTREZID", column="SYMBOL")
head(y$genes)
y <- y[!is.na(y$genes$Symbol), ]
dim(y)


#filtering low read counts
#we filter genes that have counts below 10.. cpm should be above 0.10 

keep <- rowSums(cpm(y) > 0.5) >= 2 # 2 = biological replicates 
table(keep)
y <- y[keep, , keep.lib.sizes=FALSE]


## NORMALISATION OF READS 
y<- calcNormFactors(y) #normalisation done by TMM method 
y$samples #low number of normalisation value , high no of up regulated genes. 

pch <- c(1,2,15,16,17)
colors <- rep(c("darkgreen", "red", "blue"),2 )
plotMDS(y, col = colors[group],pch = pch[group])
legend("topleft", legend=levels(group), pch=pch, col=colors, ncol=2)

plotMD(y, column=1)
abline(h=0, col="red", lty=2, lwd=2)


## DESIGN MATRIX 
design <- model.matrix(~0+group)
colnames(design) <- levels(group)
design

## DISPERSION ESTIMATION 
#edgeR uses negative binomial distribution  
y <- estimateDisp(y, design, robust = TRUE) #empirical Bayes moderated dispersion for each gene 
plotBCV(y)# Biological Coefficient of Variation 
#NB will be higher for genes with low counts. 

fit <- glmQLFit(y, design, roubust = TRUE)
head(fit$coefficients)
plotQLDisp(fit)

summary(fit$df.prior)

## DIFFERENTIAL GENE EXPRESSION 
#contrast between lactating and basal pregnant condition. 
B.Lsvp <- makeContrasts(B.lactating - B.pregnant, levels = design)
#positive LFC represents that the gene in lactating mice is upregulated relative to the pregnant mice while the negative LFC represents the pregnant mice gene is upregulated than the lactating mice 
#we use QL f-test for stricter error control 
res <- glmQLFTest(fit, contrast = B.Lsvp)
View(res)
topTags(res)


is.de <-decideTestsDGE(res) #decideTestsDGE function is used to calculate the no of up, down and not sign genes using the multiple testing and significance cut off. 
summary(is.de)


plotMD(res, status=is.de, values=c(1,-1), col=c("red","blue"),legend="topright")


tr <- glmTreat(fit, contrast=B.Lsvp, lfc=log2(1.5))
topTags(tr)


is.de <- decideTestsDGE(tr)
summary(is.de)

plotMD(tr, status=is.de, values=c(1,-1), col=c("red","blue"),legend="topright")


logCPM <- cpm(y, prior.count=2, log=TRUE)
rownames(logCPM) <- y$genes$Symbol
colnames(logCPM) <- paste(y$samples$group, 1:2, sep = "-")


o <- order(tr$table$PValue)
logCPM <- logCPM[o[1:30],]

# scale the genes to have a mean of 0 and SD of 1 
logCPM <- t(scale(t(logCPM)))

library(gplots)
col.pan <- colorpanel(100, "blue", "white", "red")
heatmap.2(logCPM, col=col.pan, Rowv=TRUE, scale="none",trace="none", dendrogram="both", cexRow=1, cexCol=1.4, density.info="none", margin=c(10,9), lhei=c(2,10), lwid=c(2,6))


##COMPARISON OF 3 OR MORE GROUPS 
#pairwise comparisons b/w groups 
con <- makeContrasts(L.PvsL = L.pregnant - L.lactating, L.VvsL = L.virgin - L.lactating, L.VvsP = L.virgin - L.pregnant, levels = design)

res <- glmQLFTest(fit, contrast = con)
topTags(res)

#test whether the change in expression between lactating and pregnant mice is the same for basal cells as it is for luminal cells
con <- makeContrasts((L.lactating-L.pregnant)-(B.lactating-B.pregnant), levels = design)
res <- glmQLFTest(fit, contrast = con)
topTags(res)



## PATHWAY ANALYSIS 
go <- goana(tr, species = "Mm")
kegg <- kegga(tr, species = "Mm")
topKEGG(kegg, n = 15, truncate = 34)


