#Group III

# alínea b)
library(survival)
library(ranger)
library(ggplot2)
library(dplyr)
library(ggfortify)
clinical_annotation <- read.csv('TCGA_COADREAD_ClinicalAnnotation.csv', header = TRUE)

clinical_annotation$status <- with(clinical_annotation,
                                   ifelse(clinical_annotation$vital_status.demographic
                                          == 'Alive', 0, 1 ))

km <- with(clinical_annotation, Surv(days_to_death.demographic, status))
plot(km)

km_cms_fit <- survfit(Surv(days_to_death.demographic, status) ~ CMS_network, data=clinical_annotation)
autoplot(km_cms_fit)

surv_diff <- survdiff(Surv(days_to_death.demographic, status) ~ CMS_network, data=clinical_annotation)
surv_diff




# Alínea C)

# devtools::install_github("Lothelab/CMScaller")
library(Biobase)
library(CMScaller)
par(mfrow=c(1,2))

data_colo <- read.csv('41598_2020_69083_MOESM2_ESM.csv')
### CMS prediction of TCGA primary colorectal cancers
res <- CMScaller(data_colo, RNAseq=TRUE, doPlot=TRUE)
head(res)

### Camera Gene Set Analysis with CMS informative gene sets
cam <- CMSgsa(emat=data_colo, class=res$prediction, RNAseq=TRUE)
head(cam$CMS4)

### limma differential gene expression analysis and visualization
deg <- subDEG(emat=data_colo, class=res$prediction, doVoom=TRUE)
subVolcano(deg, geneID="symbol")