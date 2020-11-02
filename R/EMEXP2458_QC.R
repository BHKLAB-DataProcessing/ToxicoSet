library(SummarizedExperiment)
library(ToxicoGx)
#load TSet
EMEXP2458<- readRDS("../QC/EMEXP2458.rds")

#QC1 - check correlation of replicates in each cell line
EMEXP2458_HepG2 <- ToxicoGx::subsetTo(EMEXP2458, cell_lines = "Hep-G2")

pheno_HepG2 <- as.data.frame(colData(EMEXP2458_HepG2@molecularProfiles$rna))
pheno_rep1_HepG2 <- subset(pheno_HepG2, pheno_HepG2$individual_id == "1", drop = F)
pheno_rep2_HepG2 <- subset(pheno_HepG2, pheno_HepG2$individual_id == "2", drop = F)
pheno_rep3_HepG2 <- subset(pheno_HepG2, pheno_HepG2$individual_id == "3", drop = F)

expr_mat_HepG2 <- assay(EMEXP2458_HepG2@molecularProfiles$rna)

expr_rep1_HepG2 <- expr_mat_HepG2[,colnames(expr_mat_HepG2) %in% pheno_rep1_HepG2$samplename]
expr_rep2_HepG2 <- expr_mat_HepG2[,colnames(expr_mat_HepG2) %in% pheno_rep2_HepG2$samplename]
expr_rep3_HepG2 <- expr_mat_HepG2[,colnames(expr_mat_HepG2) %in% pheno_rep3_HepG2$samplename]

cor_all_p_HepG2 <- cor(expr_rep1_HepG2, expr_rep2_HepG2, method = "pearson", use = "complete.obs")

test_cor_p_HepG2 <- cor_all_p_HepG2[order(diag(cor_all_p_HepG2), decreasing = TRUE),order(diag(cor_all_p_HepG2), decreasing = TRUE)]

png("../QC/QC1_cor_rep_diag_HepG2.png")
plot(diag(test_cor_p_HepG2), main = "Pearson correlation of HepG2 replicates", ylab = "correlation values", type = "b")
dev.off()

EMEXP2458_HepaRG <- ToxicoGx::subsetTo(EMEXP2458, cell_lines = "Hep-G2")

pheno_HepaRG <- as.data.frame(colData(EMEXP2458_HepaRG@molecularProfiles$rna))
pheno_rep1_HepaRG <- subset(pheno_HepaRG, pheno_HepaRG$individual_id == "1", drop = F)
pheno_rep2_HepaRG <- subset(pheno_HepaRG, pheno_HepaRG$individual_id == "2", drop = F)
pheno_rep3_HepaRG <- subset(pheno_HepaRG, pheno_HepaRG$individual_id == "3", drop = F)

expr_mat_HepaRG <- assay(EMEXP2458_HepaRG@molecularProfiles$rna)

expr_rep1_HepaRG <- expr_mat_HepaRG[,colnames(expr_mat_HepaRG) %in% pheno_rep1_HepaRG$samplename]
expr_rep2_HepaRG <- expr_mat_HepaRG[,colnames(expr_mat_HepaRG) %in% pheno_rep2_HepaRG$samplename]
expr_rep3_HepaRG <- expr_mat_HepaRG[,colnames(expr_mat_HepaRG) %in% pheno_rep3_HepaRG$samplename]

cor_all_p_HepaRG <- cor(expr_rep1_HepaRG, expr_rep2_HepaRG, method = "pearson", use = "complete.obs")

test_cor_p_HepaRG <- cor_all_p_HepaRG[order(diag(cor_all_p_HepaRG), decreasing = TRUE),order(diag(cor_all_p_HepaRG), decreasing = TRUE)]

png("../QC/QC1_cor_rep_diag_HepaRG.png")
plot(diag(test_cor_p_HepaRG), main = "Pearson correlation of HepaRG replicates", ylab = "correlation values", type = "b")
dev.off()
######################################################################################################################################################

#QC2 - checking ToxicoGx functions
library(data.table)
library(ToxicoGx)
png("../QC/QC2_TGx_plotting.png")
ToxicoGx::drugGeneResponseCurve(tSet = EMEXP2458, duration = c("12", "48"), cell_lines = "HepaRG", mDataTypes = "rna", 
                                dose = c("Control","High"),features = "ENSG00000000003_at",
                                drug = "Estradiol", summarize_replicates = F, verbose = T)
dev.off()

#check the distribution of normalized gene expression values
png("../QC/QC2_norm_genes.png")
hist(assay(EMEXP2458@molecularProfiles$rna), xlab = "Normalized gene expression values", main = "Distribution of normalized gene expression values")
dev.off()

######################################################################################################################################################
#QC3 - From paper
#PCA of all experimental conditions on the two cell lines. (Fig 2)
#extracting data from tset
library(affy)
library(ToxicoGx)
library(SummarizedExperiment)
library(rgl)

runENSGversions <- function(TSet){
  
  
  se <- TSet@molecularProfiles$rna
  
  expr <- assay(se)
  
  ourData <- expr
  
  geneInfo <- featureInfo(TSet, "rna")
  
  ourDataEnsgMap <- geneInfo[rownames(ourData), "gene_id"]
  
  ourDataKeep <- !is.na(ourDataEnsgMap) & !(ourDataEnsgMap %in% ourDataEnsgMap[duplicated(ourDataEnsgMap)])
  print("genes cleaned up from TSet")
  
  
  ourData <- ourData[ourDataKeep,]
  
  rownames(ourData) <- ourDataEnsgMap[ourDataKeep]
  
  colnames(ourData) <- se[,colnames(ourData)]$Array.Data.File
  
  return(ourData)
}


op_ensgversion <- runENSGversions(TSet = EMEXP2458)

se <- EMEXP2458@molecularProfiles$rna
pca.EMEXP2458 <- prcomp(t(op_ensgversion), center=TRUE, scale.=TRUE)
cell <- colData(se)[match(rownames(pca.EMEXP2458$x), se$Array.Data.File), "cellid"]
time <- colData(se)[match(rownames(pca.EMEXP2458$x), se$Array.Data.File), "duration"]

plot3d(pca.EMEXP2458$x[,1:3], col=c("HepaRG"="red", "Hep-G2" = "blue")[cell], size=6, main="By cell lines")
legend3d("topright", legend = c("HepaRG","Hep-G2") , pch = 16, col = c("HepaRG"="red", "Hep-G2" = "blue"), cex=2, inset=c(0.02))

snapshot3d(filename = '../QC/QC3_PCA_CLs.png', fmt = 'png')

plot3d(pca.EMEXP2458$x[,1:3], col=c("12"="orange", "48" = "turquoise")[as.character(time)], size=6, main="By time points")
legend3d("topright", legend = c("12hr","48hr") , pch = 16, col = c("12"="orange", "48" = "turquoise"), cex=2, inset=c(0.02))

snapshot3d(filename = '../QC/QC3_PCA_T.png', fmt = 'png')


