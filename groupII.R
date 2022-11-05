source("../plot_templates.R", local = knitr::knit_global())

library(data.table)
library(ggplot2)  
#BiocManager::install("edgeR")
genereads <- as.data.frame(colSums(read_counts[,-1]))
names(genereads)[names(genereads) == "colSums(read_counts[, -1])"] <- "Total_counts"
genereads$names <- rownames(genereads)
genereads$names <- factor(genereads$names, levels = genereads$names)
genereads <- apply(genereads,2,sort,decreasing=T)

density <- density(genereads$Total_count)
plot(density,  main = 'Density of read coverage', xlab='read counts',
     ylab='density')

#divide_by_maximum <- sweep(read_counts[,-1], 2, apply(read_counts[,-1], 2, max), "/")
#heatmap(as.matrix(divide_by_maximum[-1, ]))
set.seed(0)
lib_matrix <- read_counts[sample(nrow(read_counts[, -1]), 3925),]

#sorting columns by ascending order
lib_matrix_sorted <- apply(lib_matrix,2,sort,decreasing=T)

plotcomplexity_3 <- complexityplot_byexprgenes(data=lib_matrix_sorted,samples_to_plot=colnames(lib_matrix_sorted), 
                                               reads_in_percentage=TRUE)
plotcomplexity_3

# part B
library(edgeR)
library(limma)
library(ggplot2)

## Trying to verify the normalization

pseudo_counts <- log2(as.matrix(read_counts[c('TCGA.WS.AB45.01A', 'TCGA.T9.A92H.01A', 'TCGA.SS.A7HO.01A', 'TCGA.4T.AA8H.01A', 'TCGA.4N.A93T.01A', 'TCGA.3L.AA1B.01A' )]) + 1)
df_raw <- melt(pseudo_counts, id = rownames(read_counts))
names(df_raw)[1:2]<- c("id", "sample")
df_raw$method <- rep("Raw counts", nrow(df_raw))  

pseudo_TMM <- log2(cpm(after_normalize) + 1)

df_TMM <- melt(pseudo_TMM[, c('TCGA.WS.AB45.01A', 'TCGA.T9.A92H.01A', 'TCGA.SS.A7HO.01A', 'TCGA.4T.AA8H.01A', 'TCGA.4N.A93T.01A', 'TCGA.3L.AA1B.01A')], id = rownames(raw_counts_wn))
names(df_TMM)[1:2] <- c ("id", "sample")
df_TMM$method <- rep("TMM", nrow(df_TMM))


df_allnorm <- rbind(df_raw, df_TMM)
df_allnorm$method <- factor(df_allnorm$method,
                            levels = c("Raw counts","TMM"))

p <- ggplot(data=df_allnorm, aes(x=sample, y=value, fill=method))
p <- p + geom_boxplot()  
p <- p + theme_bw()
p <- p + ggtitle("Boxplots of normalized pseudo counts for the given samples ")
p <- p + facet_grid(. ~ method) 
p <- p + ylab(expression(log[2] ~ (normalized ~ count + 1))) + xlab("")
p <- p + theme(title = element_text(size=10), axis.text.x = element_blank(), 
               axis.ticks.x = element_blank())
print(p)

##
#colData <- data.frame(condition=ifelse(aa == 0, "control", "triggered"))
#dds <- DESeqDataSetFromMatrix(round(ReadCounts), colData = colData, formula(~ condition))
#dds <- estimateSizeFactors(dds)
#normalized_counts <- counts(dds, normalized=TRUE)


dge <- DGEList(read_counts[,-1])
dge <- calcNormFactors(dge, method = "TMM")
v <- voom(dge,plot=FALSE)

# part C
library(stringr)
library(tibble)
library(patchwork)
pca =  prcomp(t(v$E), center = TRUE,scale. = TRUE)
screeplot(pca, main="PCA Analysis", xlab="PCi")

pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,1],
                       Y=pca$x[,2])

pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100,2)

dev.new()
ggplot(pca.data, aes(x=X, y=Y)) +
  geom_point() +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  theme_bw() +
  ggtitle('PCA Analysis')

#dev.new()
label1 <-cut(genereads$Total_counts, c(0, 1*10^7, 3*10^7, 5*10^7, 7*10^7, 1*10^8))

pca_data <- data.frame(pca$x)
dataplot1 <- data.frame(label1, pca_data[,1], pca_data[,2])
colnames(dataplot1)<- c("label","x","y")
plot1 <- ggplot(dataplot1, aes(x = x, y = y, color = label)) + 
  geom_point()+
  guides(color=guide_legend("Total No. of Reads")) +
  xlab(label="PC1") + 
  ylab(label = "PC2") +
  ggtitle("PCA Analysis per Read Coverage")



### Clinical study
clinical_info <- read.delim(file = '../TCGA_COADREAD_ClinicalAnnotation.txt', 
                            header=TRUE,stringsAsFactors = FALSE)

colnames_lib <- colnames(read_counts[-1])
colnames_t <- str_sub(colnames(read_counts[-1]), end=12) 
colnames_t<- gsub('\\.', '-', colnames_t)
colnames_index<-match(colnames_t,str_sub(clinical_info[,1], ))

new_clinical = clinical_info[FALSE,]

i_set<- seq(1, length(colnames_t))
for (i in i_set) {
  new_clinical[i,]<-clinical_info[colnames_index[i],]
}

loading_scores <-pca$rotation[,1]
gene_scores <- abs(loading_scores)
gene_score_ranked <- sort(gene_scores, decreasing=TRUE)
top_genes <- names(gene_score_ranked[1:20])

top_genes 

pca$rotation[top_genes,1]

loading_scores <-pca$rotation[,2]
gene_scores <- abs(loading_scores)
gene_score_ranked <- sort(gene_scores, decreasing=TRUE)
top_genes <- names(gene_score_ranked[1:20])

top_genes 

pca$rotation[top_genes,2]

# Genero

label1 <-new_clinical$gender.demographic

dataplot1 <- data.frame(label1, pca_data[,1], pca_data[,2])
colnames(dataplot1)<- c("label","x","y")

#dev.new()
plot2 <- ggplot(dataplot1, aes(x = x, y = y, color = label)) + 
  geom_point()+
  guides(color=guide_legend("Gender")) +
  xlab(label="PC1") + 
  ylab(label = "PC2") +
  ggtitle("PCA Analysis per Gender")

# Idade do diagnostico

label1 <-new_clinical$age_at_initial_pathologic_diagnosis

dataplot1 <- data.frame(label1, pca_data[,1], pca_data[,2])
colnames(dataplot1)<- c("label","x","y")

#dev.new()
plot3<- ggplot(dataplot1, aes(x = x, y = y, color = label)) + 
  geom_point()+
  guides(color=guide_legend("Age at Diagnosis")) +
  xlab(label="PC1") + 
  ylab(label = "PC2") +
  ggtitle("PCA Analysis per Age at Diagnosis")


# Type of sample 

label1 <-new_clinical$aa

dataplot1 <- data.frame(label1, pca_data[,1], pca_data[,2])
colnames(dataplot1)<- c("label","x","y")

#dev.new()
plot4 <- ggplot(dataplot1, aes(x = x, y = y, color = label)) + 
  geom_point()+
  guides(color=guide_legend("Type of Sample")) +
  scale_color_manual(labels = c("Healthy","Primary", "Metastasis"),values = c ("purple","pink","black"))+
  xlab(label="PC1") + 
  ylab(label = "PC2") +
  ggtitle("PCA Analysis per Type of Sample")

# Part D
library(ggrepel)
library(clusterProfiler)
library(enrichplot)

aa <- str_sub(colnames_lib, start=14, end = 15)
for (i in i_set){
  if (aa[i]=="01" | aa[i] == "02") #tumor
  {aa[i]<- as.numeric(1) }
  else if (aa[i]=="06") #Metastatic
  {aa[i]<- as.numeric(2)}
  else {aa[i]<- as.numeric(0)} #healthy
}

new_clinical <- add_column(new_clinical, aa, .after = "patientID")
new_clinical <- add_column(new_clinical, genereads$Total_counts, .after = "aa")

design.matrix <- model.matrix(~aa*age_at_initial_pathologic_diagnosis, data = new_clinical)

read_counts <- read_counts[rowSums(read_counts[, -1]) > 0, ]
rownames(read_counts)<- read_counts[,1]
read_counts<- read_counts[,-1]

read_counts <- subset(read_counts, select = -c(TCGA.5M.AATA.01A, 
                                               TCGA.5M.AAT5.01A, 
                                               TCGA.F5.6810.01A))

dge2 <- DGEList(counts=data.matrix(read_counts))
dge2 <- calcNormFactors(dge2)
v2 <- voom(dge2,design.matrix,plot=TRUE)
linearfit = lmFit(v2$E,design.matrix) 
eBfit = eBayes(linearfit)
volcanoplot(eBfit,coef=2,style="B-statistic")
mypval=0.05
myfc=2

#tumorvsnormal
limma.res.pval <- topTable(eBfit,coef=2,sort.by = "t", number = Inf, p.val=mypval)
dim(limma.res.pval)

limma.res.pval.FC <- limma.res.pval[which(abs(limma.res.pval$logFC)>myfc),]
dim(limma.res.pval.FC)

#normalvstumorvsage
limmaage.res.pval <- topTable(eBfit,coef=5, sort.by = "t", number = Inf, p.val=mypval)
dim(limmaage.res.pval)

limmaage.res.pval.FC <- limmaage.res.pval[which(abs(limmaage.res.pval$logFC)>0.05),]
dim(limmaage.res.pval.FC)

HealthyTumor <- topTable(eBfit,coef=2,sort.by = "t", number = Inf)
HealthyTumorAge <- topTable(eBfit,coef=5,sort.by = "t", number = Inf)

##############
#Healthy and Tumor
#############

#now for plotting --- a fancier version of this
#volcanoplot(eBfit,coef=2,style="B-statistic")

#creating a column with p values, where we deem significant if p-value under 5%
significance_ht <- HealthyTumor$adj.P.Val

i_set = seq(1, length(significance_ht))
for (i in i_set){
  if (as.numeric(significance_ht[i]) < 0.05){
    significance_ht[i] <- "Significant"
  } else {
    significance_ht[i] <- "Not significant"
  }
}
significance_ht<-data.frame(significance_ht)

#create data frame for easier plotting
plot_ht <- data.frame(significance_ht, HealthyTumor$logFC, HealthyTumor$B)
colnames(plot_ht)<- c("label","x","y")
#add column with gene names
plot_ht <- cbind(gene=rownames(HealthyTumor), plot_ht) 
limma.FC <- cbind(gene=rownames(limma.res.pval.FC), limma.res.pval.FC) 
interm_fc<-(match(plot_ht$gene,limma.FC$gene, nomatch = NA))
gene_labels <- data.frame(plot_ht$gene)
plot_ht <- cbind(glabels = gene_labels, plot_ht) 

dev.new()
ggplot(plot_ht, aes(x = x, y = y,color = label)) + #volcanoplot with log2Foldchange versus B statistic
  geom_point() + #add points colored by significance
  guides(color=guide_legend("Significance")) +
  xlab(label="Log2FoldChange") + 
  ylab(label = "Log Odds of Differential Expression") +
  ggtitle("Difference between healthy and tumor") +
  geom_text_repel(data=head(plot_ht, 30), aes(label=plot_ht.gene)) #plot 20 most diff. expressed genes
#scale_color_manual(values=c("black", "blue")) 

################
#Healthy and Tumor and Age
################

#creating a column with p values, where we deem significant if p-value under 5%
significance_ht <- HealthyTumorAge$adj.P.Val

i_set = seq(1, length(significance_ht))
for (i in i_set){
  if (as.numeric(significance_ht[i]) < 0.05){
    significance_ht[i] <- "Significant"
  } else {
    significance_ht[i] <- "Not significant"
  }
}
significance_ht<-data.frame(significance_ht)

#create data frame for easier plotting
plot_ht <- data.frame(significance_ht, HealthyTumorAge$logFC, HealthyTumorAge$B)
colnames(plot_ht)<- c("label","x","y")
#add column with gene names
plot_ht <- cbind(gene=rownames(HealthyTumorAge), plot_ht) 
limma.FC <- cbind(gene=rownames(limmaage.res.pval.FC), limmaage.res.pval.FC) 
interm_fc<-(match(plot_ht$gene,limma.FC$gene, nomatch = NA))
gene_labels <- data.frame(limma.FC$gene[interm_fc])
plot_ht <- cbind(glabels = gene_labels, plot_ht) 

dev.new()
ggplot(plot_ht, aes(x = x, y = y,color = label)) + #volcanoplot with log2Foldchange versus B statistic
  geom_point() + #add points colored by significance
  guides(color=guide_legend("Significance")) +
  xlab(label="Log2FoldChange") + 
  ylab(label = "Log Odds of Differential Expression") +
  ggtitle("Healthy, Tumor and Age ") +
  geom_text_repel(data=head(plot_ht, 98), aes(label=limma.FC.gene.interm_fc.)) #plot 20 most diff. expressed genes
#scale_color_manual(values=c("black", "blue")) 


### pathway analyse
choose_analyse <- HealthyTumorAge

original_gene_list <- choose_analyse$t

# name the vector
names(original_gene_list) <- rownames(choose_analyse)

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

# using KEGG
ids<-bitr(names(original_gene_list), 
          fromType = "SYMBOL", toType = "ENTREZID", OrgDb=organism)

dedup_ids = ids[!duplicated(ids[c("SYMBOL")]),]

df2 = choose_analyse[rownames(choose_analyse) %in% dedup_ids$SYMBOL,]

df2$Y = dedup_ids$ENTREZID
kegg_gene_list <- df2$t
names(kegg_gene_list) <- df2$Y

kegg_gene_list<-na.omit(kegg_gene_list)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)

kegg_organism = "hsa"

kk2 <- gseKEGG(geneList     = kegg_gene_list, organism = kegg_organism,
               nPerm        = 10000, minGSSize    = 3, maxGSSize    = 800,
               pvalueCutoff = 0.05, pAdjustMethod = "none", 
               keyType       = "ncbi-geneid")
dev.new()
dotplot(kk2, showCategory = 10, title = "Enriched Pathways - Healthy vs Tumor vs Age" , split=".sign") + facet_grid(.~.sign)
