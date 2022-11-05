#BiocManager::install("biomaRt")
library(biomaRt)
library(dplyr)

abundance <- read.table('abundance.tsv', header=TRUE)

read_counts <- read.table('../TCGA_COADREAD_Gene_ReadCounts.txt', header = TRUE)
read_counts <-read_counts[,c("Gene","TCGA.DM.A288.01A")]
  
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

transcript_ids <-abundance$target_id

res <- getBM(attributes = c('ensembl_transcript_id_version', 
                            'external_gene_name'),
             filters = 'ensembl_transcript_id_version', 
             values = transcript_ids,
             mart = mart)


data <- merge(res, abundance, by.x="ensembl_transcript_id_version", 
              by.y = "target_id")
data <- merge(data, read_counts, by.x='external_gene_name', by.y = 'Gene')
data <- data[,c('external_gene_name', 'ensembl_transcript_id_version',
                'est_counts', 'TCGA.DM.A288.01A')]


summary = data %>% group_by(external_gene_name)  %>% 
            summarise(est_counts = sum(est_counts), .groups = 'drop')

comparing <- merge(summary, read_counts,  by.x='external_gene_name', by.y = 'Gene')

plot(comparing$est_counts, comparing$TCGA.DM.A288.01A,  xlab = "kallisto counting", 
     ylab = "TCGA-DM-A288-01A counting", 
     main = 'Comparing between two methods',  col = "blue", log='xy',
     xlim=range(1:5e+5) 
    )
abline(b=1, a=-0.1, col="red")
legend(1, 1.5e+5, legend="a=-0.1 ; b=1",
       col=c("red", "blue"), lty=1:2, cex=0.8)
#abline(lm(comparing$est_counts ~ comparing$TCGA.DM.A288.01A))