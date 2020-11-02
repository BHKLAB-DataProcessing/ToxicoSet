
## QC FOR TGH - Plotting gene expression changes from TG-GATEs Human dataset
  
#   Our first QC is based on a plot from a paper by Rangel-Escareño et al., in which gene expression changes of CYP1A1 (gene associated with xenobiotic metabolism) has been plotted at all concentrations and time points. The plot shows clear differential expression at time 8(hr) suggesting that higher the dose, larger the impact of CCL4 on this gene.
# For plotting the gene expression under same conditions using the package, the first step is to load the datasets from disk or download them using the downloadTSet function above. In the following example, we use the toy dataset provided with the package to illustrate the process. 
# To plot, the function drugGeneResponseCurve has been used wherein mandatory inputs such as dataset, drug name, cell-line, molecular type, gene name, dose and time points should be specified.

#**Study 1 - Plot time dependent dose response of Carbon tetra chloride on CYP1A1 gene.**
  
library(ToxicoGx)
library(ggplot2)
# Load the tset 
TGGATES_humanldh <- readRDS("../results/TGGATES_humanldh.rds")

png("../results/QC1_TGH.png")
drugGeneResponseCurve(tSet = TGGATES_humanldh, duration = c("2", "8", "24"), 
                      cell_lines = "Hepatocyte", mDataTypes = "rna", 
                      features = "ENSG00000140465_at", 
                      dose = c("Control", "Low", "Middle","High"),
                      drug = "Carbon tetrachloride", 
                      summarize_replicates = F,
                      ggplot_args=list(scale_color_manual(values=c("red", "green", "dark blue", "turquoise"))),
                      verbose = T)
dev.off()
#**Study 2 - PMID - 25399406: Fig 1A - PCA was performed using the 100 top-ranking genes with the highest fold change (absolute values) across all compounds. 2 main clusters—the lower cluster subdivided into several sub-clusters and the majority of the treated samples that move in the direction of the first principal component. .**
library(ToxicoGx)
library(car)
library(SummarizedExperiment)
#Genes of interest
probes <- readRDS("../data/QC/PCA_probes.rds")

#extracting data from tset
#featureData
feat_data <- as.data.frame(rowData(TGGATES_humanldh@molecularProfiles$rna))
#metaData
pheno_data <- as.data.frame(colData(TGGATES_humanldh@molecularProfiles$rna))
#expression values
assay_data <- assay(TGGATES_humanldh@molecularProfiles$rna)
rownames(assay_data) <- gsub("_at", "", rownames(assay_data))

#subsetting the samples - Control, high dose at 24 hr time point
samples <- subset(pheno_data$samplename, (pheno_data$dose_level == "Control" | pheno_data$dose_level == "High") & pheno_data$duration == "24")

control <- subset(pheno_data,pheno_data$dose_level=="Control" & pheno_data$duration == "24",select=c(samplename))

high <- subset(pheno_data,pheno_data$dose_level=="High" & pheno_data$duration == "24",select=c(samplename))

expr <- as.data.frame(t(subset(assay_data,rownames(assay_data) %in% probes,select=as.character(samples))))


expr$control<-NA
expr$control[rownames(expr) %in% as.character(control$samplename)] <- 1 #Control
expr$control[rownames(expr) %in% as.character(high$samplename)] <- 2 #High

#PCA & plotting
tset.pca <- prcomp(as.matrix(expr),scale. = TRUE)

png("../results/QC2_TGH.png")

scatterplot(x = tset.pca$x[,1], 
            y = tset.pca$x[,2], 
            regLine = FALSE, 
            smooth = FALSE, 
            boxplots = FALSE, 
            groups = expr$control, 
            col = c('dark green','chartreuse1'),
            cex = 1, 
            pch = c(20,20,20),
            legend=FALSE,
            xlab="PC1",
            ylab="PC2",
            main="Principal Component Analysis : High Dose, 24h")
legend(10,5,legend=c("Control","treated"),col=c('dark green','light green'),pch=c(20,20,20),cex=1, bty = 0, pt.cex = 2)
dev.off()

#**Study 3 - PMID - 30426165 : fig 4 - Effect of azathioprine on the NRF2-associated gene module(#325)**
  library(ToxicoGx)
  library(SummarizedExperiment)
  
  # Load the tset 
  drug <- "Azathioprine"
  conc <- 72.8
  
  #Genes of interest
  genes <-c("ENSG00000159231_at","ENSG00000125037_at", "ENSG00000102393_at","ENSG00000181019_at", "ENSG00000109854_at", "ENSG00000164220_at")
  
  #apply function to extract exprs matrix
  values_se <- lapply(genes, function(gene){
    #subset pehnodata for desired drugs
    drug_subset <- subset(as.data.frame(colData(TGGATES_humanldh@molecularProfiles$rna)),drugid == drug,select=c(samplename, dose_level, individual_id,concentration))
    #subset for ony high conc
    drug_subset_high <- subset(drug_subset, concentration == conc)
    #extracting exprs
    assay <- assay(TGGATES_humanldh@molecularProfiles$rna)
    #subsetting exprs matrix
    drug_subset$expression <- assay[gene,as.character(drug_subset$samplename)]
    drug_subset_high$expression <- assay[gene,as.character(drug_subset_high$samplename)]
    #ctrl rep
    ctrlA <- na.omit(drug_subset$expression[select=c(drug_subset$dose_level == "Control" & drug_subset$individual_id=="1")])
    ctrlB <- na.omit(drug_subset$expression[select=c(drug_subset$dose_level == "Control" & drug_subset$individual_id=="2")])
    
    highA <- na.omit(drug_subset_high$expression[select=c(drug_subset_high$dose_level == "High" & drug_subset_high$individual_id=="1")])
    highB <- na.omit(drug_subset_high$expression[select=c(drug_subset_high$dose_level == "High" & drug_subset_high$individual_id=="2")])
    
    ctrl <- rowMeans(cbind(ctrlA, ctrlB))
    high <- rowMeans(cbind(highA, highB))
    
    
    normalised_vehicle <- (high-ctrl)*100
    return(normalised_vehicle)
  })
  
  values_se <- as.data.frame(do.call(rbind,values_se))
  colnames(values_se) <- c(2,8,24)
  rownames(values_se) <- genes
  
  
  time <- c(2,8,24)
  legendnames <- c('CBR3','EMC3','GLA','NQO1','HTAT1P2','F2RL2')
  colours <- c("purple", "red","green","violet","orange", "turquoise")
  
  png("../results/QC3_TGH.png")
  
  matplot(x = time, y = t(values_se)+100, col=colours, 
          pch=rep(21,ncol(values_se)), type=c("b"), lty=rep(1,ncol(values_se)), lwd=rep(5,ncol(values_se)),
          bg=colours,
          xlim=range(0,4,8,12,16,20,24),ylim=range(0,100,200,300,400,500),main="Azathioprine (Mod 325)",
          xlab="Time",ylab="mRNA Level (% Vehicle)")
  par(font=2) 
  legend(23.1,520,legend=legendnames, 
         col = colours,
         pch=rep(21,ncol(values_se)), pt.bg = colours,
         text.col = 'black',
         lty=rep(1,ncol(values_se)),lwd = rep(5, ncol(values_se)), cex=0.75,xjust = 0.5,bty = "n", adj = 0.25)
  
  dev.off()
  
  ## Connectivity map analysis on TG-GATEs and human hepatocarcinoma signatures.
  
  # For the second case study, we will recreate an analysis from the paper by Jos Kleinjans et al., wherein connectivity mapping has been used to predict compound carcinogenicity
  # by linking in vivo human hepatocarcinoma (HCC) signature geneset with in vitro TG-GATEs primary human hepatocyte data. In this example, we are using the toy dataset. The full dataset has to be downloaded to carry out the whole analysis done in the paper.
  # The HCC signature, already mapped to the gene level, has been included in this package and it can be loaded by calling data(HCC_sig). Once the dataset is loaded, recreate drug signatures for each drug using the function drugPerturbationSig to perform statistical modelling of the transcriptomic response to the application of each drug. We then compare the observed up-regulated and down-regulated genes to HCC signature published in the paper. The output will be the GSEA connectivity score with FDR values that can be used to determine the correlation between the two signatures.
  
#  **Study 4 - Case study - PMID: 23940306 **
    
  require(xtable)
  library(ToxicoGx)
  library(xtable)
  library(Biobase)
  library(car)
  library(ggplot2)
  library(PharmacoGx)
  library(SummarizedExperiment)
  
  
  #HCC signature genes downloaded from supp data.
  raw_hcc <- read.delim("../data/QC/Supp._data_2.txt", sep = "\t") 
  colnames(raw_hcc) <- c("feature","direction")
  
  ff <- as.data.frame(rowData(TGGATES_humanldh@molecularProfiles$rna))
  
  #merge fdata with raw hcc to map entrez id
  
  merge_entrez <- merge(raw_hcc, ff, by.x = "feature", by.y = "EntrezGene.ID", all.x = T)
  
  ss <- subset(merge_entrez, select = c("gene_id","direction"), drop= F)
  
  colnames(ss)[1] <- "feature"
  HCC_sig <- subset(ss, subset = !is.na(ss$feature), drop = F)
  
  HCC_sig <- HCC_sig[order(HCC_sig$direction),]
  
  rownames(HCC_sig) <- HCC_sig[,1]
  
  
  
  # drug.perturbation <- ToxicoGx::drugPerturbationSig(tSet = TGGATES_humanldh,mDataType="rna",cell_lines = "Hepatocyte",
  #                                           duration = "24", dose = c("Control", "Low"),drugs = ToxicoGx::drugNames(TGGATES_humanldh),verbose=FALSE)
  #                                          
  #saveRDS(drug.perturbation, "../data/QC/drug.perturbationAllCtrlLow.rds")
  
  drug.perturbation <- readRDS("../data/QC/drug.perturbationAllCtrlLow.rds")
  
  # res_all <- apply(drug.perturbation[,,c("tstat", "fdr")],
  #                  2, function(x, HCC){
  #                    return(connectivityScore(x=x,
  #                                             y=HCC[,2,drop=FALSE],
  #                                             method="fgsea", nperm=1000, nthread = 4))
  #                  }, HCC=HCC_sig)
  # 
  # rownames(res_all) <- c("Connectivity", "P Value")
  # res_all <- t(res_all)
  # res_all <-  cbind(res_all,"FDR" = p.adjust(res_all[,2], method="fdr"))
  # res_all <- res_all[order(res_all[,3]),]
  # 
  # View(res_all)
  #saveRDS(res_all, "../data/QC/res_all.rds")
  
  res_all <- readRDS("../data/QC/res_all.rds")
  xtable(res_all,
         caption='Connectivity Score results for HCC and TG-GATEs PHH gene signature.')
  
  #checking the correlation between published and recomputed scores
  sscmap.scores <- read.csv("../data/QC/CS3_table3.csv", stringsAsFactors = F)
  
  badchars <- "[\xb5]|[]|[ ,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]|[(]|[)]"
  
  sscmap.scores$cleannames <- toupper(gsub(badchars, "", sscmap.scores$Compound_changed))
  
  
  res_new <- cbind(res_all, "cleannames" = toupper(gsub(badchars, "", rownames(res_all))))
  
  mer_ssc_res <- merge(sscmap.scores, res_new, by = "cleannames", all.x = T )  
  
  mer_ssc_res <- subset(mer_ssc_res, subset = !is.na(mer_ssc_res$Connectivity), drop = F)
  
  
  sub_mer_ssc_res <- subset(mer_ssc_res, select = c("Compound_changed","setscore", "Connectivity"))
  
  sub_mer_ssc_res$Connectivity <- as.numeric(as.character(sub_mer_ssc_res$Connectivity))
  
  #check spearman's correlation
  cor.test(sub_mer_ssc_res$Connectivity,sub_mer_ssc_res$setscore, method = "s")
  
  #plot the correlation
  #pdf("CS3_corrplot_all.pdf", width = 25, height = 20)
  png("../results/QC4.1_TGH.png")
  
  scatterplot(sub_mer_ssc_res$Connectivity,sub_mer_ssc_res$setscore,boxplots = F, grid = T, xlab = "ConnectivityScore", ylab = "sscMAPscore",cex = 1,
              pch = 20, cex.lab=1, cex.axis=1, cex.main=1.5, cex.sub=1.5, par(mar = c(7, 7, 0, 0)))
  
  dev.off()
  
  #filter out FDR > 0.1
  res_all_fil <- subset(res_all, subset = res_all[,"FDR"] <  0.1, drop =F)
  
  badchars <- "[\xb5]|[]|[ ,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]|[(]|[)]"
  
  sscmap.scores$cleannames <- toupper(gsub(badchars, "", sscmap.scores$Compound_changed))
  
  
  res_new_fil <- cbind(res_all_fil, "cleannames" = toupper(gsub(badchars, "", rownames(res_all_fil))))
  
  mer_ssc_res_fil <- merge(sscmap.scores, res_new_fil, by = "cleannames", all.x = T )  
  
  mer_ssc_res_fil <- subset(mer_ssc_res_fil, subset = !is.na(mer_ssc_res_fil$Connectivity), drop = F)
  
  
  sub_mer_ssc_res_fil <- subset(mer_ssc_res_fil, select = c("Compound_changed","setscore", "Connectivity"))
  
  sub_mer_ssc_res_fil$Connectivity <- as.numeric(as.character(sub_mer_ssc_res_fil$Connectivity))
  
  #check spearman's correlation
  cor.test(sub_mer_ssc_res_fil$Connectivity,sub_mer_ssc_res_fil$setscore, method = "s")
  
  #plot the correlation
  #pdf("CS3_corrplot_fdrfilter.pdf", width = 20, height = 15)
  png("../results/QC4.2_TGH.png")
  
  scatterplot(sub_mer_ssc_res_fil$Connectivity,sub_mer_ssc_res_fil$setscore,boxplots = F, grid = T, xlab = "ConnectivityScore", ylab = "sscMAPscore",cex = 1,
              pch = 20, cex.lab=1, cex.axis=1, cex.main=1.5, cex.sub=1.5)
  dev.off()
  
#  In the table, certain drugs show a positive connectivity score. This observation aligns with the trends reported in the paper. 
  
  
  
  ## QC FOR TGR - Comparing with old TSets (No published studies used for QC due to lack of suitability)
  
  
#  **Study 1 - Plot time dependent dose response of Carbon tetra chloride on Cyp1a1 gene.**
  library(ToxicoGx)
  # Load the tset 
  TGGATES_ratldh <- readRDS("../results/TGGATES_ratldh.rds")
  png("../results/QC1_TGR.png")
  
  drugGeneResponseCurve(tSet = TGGATES_ratldh, duration = c("2", "8", "24"), 
                        cell_lines = "Hepatocyte", mDataTypes = "rna", 
                        features = "ENSRNOG00000019500_at", 
                        dose = c("Control", "Low", "Middle","High"),
                        drug = "Carbon tetrachloride", 
                        summarize_replicates = F,
                        ggplot_args=list(scale_color_manual(values=c("green", "blue", "violet", "red"))),
                        verbose = T)
  
  dev.off()

  
#  **Study 2 - rat counterpart study of TGH Study 3.**
  library(ToxicoGx)
  library(SummarizedExperiment)
  
  # Load the tset 
  drug <- "Azathioprine"
  conc <- 3.6
  
  #Genes of interest
  genes <-c("ENSRNOG00000001701_at","ENSRNOG00000009934_at", "ENSRNOG00000012772_at","ENSRNOG00000022189_at", "ENSRNOG00000018054_at")
  
  #apply function to extract exprs matrix
  values_se <- lapply(genes, function(gene){
    #subset pehnodata for desired drugs
    drug_subset <- subset(as.data.frame(colData(TGGATES_ratldh@molecularProfiles$rna)),drugid == drug,select=c(samplename, dose_level, individual_id,concentration))
    #subset for ony high conc
    drug_subset_high <- subset(drug_subset, concentration == conc)
    #extracting exprs
    assay <- assay(TGGATES_ratldh@molecularProfiles$rna)
    #subsetting exprs matrix
    drug_subset$expression <- assay[gene,as.character(drug_subset$samplename)]
    drug_subset_high$expression <- assay[gene,as.character(drug_subset_high$samplename)]
    #ctrl rep
    ctrlA <- na.omit(drug_subset$expression[select=c(drug_subset$dose_level == "Control" & drug_subset$individual_id=="1")])
    ctrlB <- na.omit(drug_subset$expression[select=c(drug_subset$dose_level == "Control" & drug_subset$individual_id=="2")])
    
    highA <- na.omit(drug_subset_high$expression[select=c(drug_subset_high$dose_level == "High" & drug_subset_high$individual_id=="1")])
    highB <- na.omit(drug_subset_high$expression[select=c(drug_subset_high$dose_level == "High" & drug_subset_high$individual_id=="2")])
    
    ctrl <- rowMeans(cbind(ctrlA, ctrlB))
    high <- rowMeans(cbind(highA, highB))
    
    
    normalised_vehicle <- (high-ctrl)*100
    return(normalised_vehicle)
  })
  
  values_se <- as.data.frame(do.call(rbind,values_se))
  colnames(values_se) <- c(2,8,24)
  rownames(values_se) <- genes
  
  
  time <- c(2,8,24)
  legendnames <- c('CBR3','EMC3','NQO1','HTAT1P2','F2RL2')
  colours <- c("purple", "red","green","violet","orange", "turquoise")
  
  png("../results/QC2_TGR.png")
  
  matplot(x = time, y = t(values_se)+100, col=colours, 
          pch=rep(21,ncol(values_se)), type=c("b"), lty=rep(1,ncol(values_se)), lwd=rep(5,ncol(values_se)),
          bg=colours,
          xlim=range(0,4,8,12,16,20,24),ylim=range(0,100,200,300,400,500),main="Azathioprine (Mod 325)",
          xlab="Time",ylab="mRNA Level (% Vehicle)")
  par(font=2) 
  legend(23.1,520,legend=legendnames, 
         col = colours,
         pch=rep(21,ncol(values_se)), pt.bg = colours,
         text.col = 'black',
         lty=rep(1,ncol(values_se)),lwd = rep(5, ncol(values_se)), cex=0.75,xjust = 0.5,bty = "n", adj = 0.25)
  
  dev.off()

  ## QC FOR DM - Comparing with old TSets (No published studies used for QC due to lack of suitability)
  
  
#  **Study 1 - Plot time dependent dose response of Carbon tetra chloride on Cyp1a1 gene.**
  library(ToxicoGx)
  # Load the tset 
  drugMatrix <- readRDS("../results/drugMatrix.rds")
  png("../results/QC1_DM.png")
  
  drugGeneResponseCurve(tSet = drugMatrix, duration = c("16", "24"), 
                        cell_lines = "Hepatocyte", mDataTypes = "rna", 
                        features = "ENSRNOG00000019500_at", 
                        dose = c("Control","High"),
                        drug = "Carbon tetrachloride", 
                        summarize_replicates = F,
                        ggplot_args=list(scale_color_manual(values=c("green", "blue", "violet", "red"))),
                        
                        verbose = T)
  
 dev.off()
  
  
  
  
  
  