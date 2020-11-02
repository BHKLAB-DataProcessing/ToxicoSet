#devtools::install_github("bhklab/CoreGx", ref = "master", force = T)
#devtools::install_github("bhklab/ToxicoGx", ref = "master", force = T)
library(ToxicoGx)
library(Biobase)
library(affy)
library(affyio)
library(BiocManager)
library(dplyr)
library(biomaRt)
library(ggplot2)
library(data.table)
library(SummarizedExperiment)

create_phenodata_DM <- function(species=c("Rat"), verbose = TRUE) { 
  if (verbose) {message("Creating phenodata_DM object...")} 
  phenodata_DM <- read.csv("../data/DM/phenodataDMSOcontrols.csv", stringsAsFactors = FALSE)
  phenodata_DM$Dose_Unit <- "μM"
  
  if(verbose) {message("phenodata_DM object created!")}
  return(phenodata_DM)
}
phenodata_DM <- create_phenodata_DM("Rat")
rownames(phenodata_DM) <- phenodata_DM$samplename
phenodata_DM$duration <- floor(phenodata_DM$duration)

#library(rat2302rnensgcdf)

create_exprsdata_DM <- function(species=c("Rat"), phenodata_DM, verbose = TRUE) {
  if (verbose) {message("Creating eset object")}
  if(species == "Rat"){
  }
  #celFiles <- paste("/DrugMatrix_raw_files/drugmatrix.rathepatocyte.celfiles","/", phenodata_DM[,"celfilename"], sep="")
  
  #esetNorm <- just.rma(filenames = celFiles, verbose = TRUE, cdfname = "rat2302rnensgcdf")
  #saveRDS(esetNorm, "data/DM/eset_DM.rds")
  eset <- readRDS("../data/DM/eset_DM.rds")
  
  storageMode(eset)<-"environment"
  eset <-subset(eset, substr(rownames(eset@assayData$exprs), 0, 4) != "AFFX")
  storageMode(eset)<-"lockedEnvironment"
  annotation(eset)<-"rna"
  if (verbose) {message("eset object created!")}
  return(eset)
}
eset <- create_exprsdata_DM("Rat", phenodata_DM)
colnames(eset) <- gsub(".CEL", "", colnames(eset))

create_featuredata_DM <- function(species=c("Rat"), eset, verbose = TRUE){
  if (verbose) {message("Creating featuredata_DM object...")}
  if (species == "Rat"){
    ensembl <- useMart("ensembl")
    datasets <- listDatasets(ensembl)
    ensembl = useMart("ensembl", dataset = "rnorvegicus_gene_ensembl", host="uswest.ensembl.org",ensemblRedirect = FALSE)
    storageMode(eset) <- "environment"
    affxrows <- rownames(eset@assayData$exprs)
    rownames(eset@assayData$exprs) <- substr(rownames(eset@assayData$exprs), 1, nchar(affxrows)-3)
    CELgenes <- affxrows
    CELgenes1 <- gsub(".at", " ", CELgenes)
    results <-getBM(attributes=c("external_gene_name","ensembl_gene_id","gene_biotype","entrezgene_id"
                                 ,"external_transcript_name","ensembl_transcript_id"), filters = "ensembl_gene_id"
                                  ,values=CELgenes1, mart=ensembl,checkFilters = TRUE)
    uniqueB <- results[!duplicated(results$ensembl_gene_id),]
    CELnotB <- unique(CELgenes1) [!unique(CELgenes1) %in% uniqueB$ensembl_gene_id]
    names(uniqueB) <- c("gene_name", "gene_id", "gene_biotype", "EntrezGene.ID", "transcript_name", "transcript_id")
    finalFeature <- uniqueB
    
    finalFeature$BEST <- NA
    names(finalFeature) <- c("Symbol", "gene_id", "gene_biotype", "EntrezGene.ID", "transcript_name", "transcript_id", "BEST")
    rownames(finalFeature) <- finalFeature$gene_id
    finalFeature$gene_id
    geneid1 <- finalFeature$gene_id
    
    for (i in 1:length(geneid1)) {
      geneid1[i] = paste(geneid1[i], "at", sep="_")
    }
    geneid1
    finalFeature$gene_id <- geneid1
    finalFeature$gene_id
    finalFeature[,1]
    rownames(finalFeature) = finalFeature$gene_id
    
    if(verbose) {message("featuredata_DM object created!")}
    return(finalFeature)
    
  }
}
featuredata_DM <- create_featuredata_DM("Rat", eset)


create_Expressionset <- function(species=c("Rat"), eset, verbose = TRUE){
  if (verbose) {message("Creating expressionset...")}
  if (species == "Rat"){

pData(eset) <- phenodata_DM
fData(eset) <- featuredata_DM
#sorting rownames to maintain feature data mapping that m=is otherwise shuffled after converting to SE
fData(eset) <- fData(eset)[sort(rownames(fData(eset))),]
stopifnot(all(rownames(fData(eset)) == rownames(exprs(eset))))
stopifnot(all(rownames(pData(eset)) == colnames(exprs(eset))))

storageMode(eset) <- "lockedEnvironment"
return(eset)
  }
}

ExpressionSet <- create_Expressionset("Rat", eset)

#the conversion function might be incorporated to the package later. This step needs to be updated then
print("Creating summarized experiment object...")
#new_SE_DM <-  as(ExpressionSet, value="SummarizedExperiment")
new_SE_DM <- SummarizedExperiment::makeSummarizedExperimentFromExpressionSet(ExpressionSet)
stopifnot(all(rownames(colData(new_SE_DM)) == rownames(pData(ExpressionSet))))
stopifnot(all(rownames(rowData(new_SE_DM)) == rownames(fData(ExpressionSet))))

print("Done!")

######################toxicoset constructor function ######################
create_curationDrug <- function(phenodata_DMData, verbose = TRUE){
  curationDrug <- unique(subset(phenodata_DM, select=c(drugid, dataset_drugid)))
  rownames(curationDrug) <- curationDrug$drugid
  names(curationDrug) <- c("unique.drugid", "dataset_drugid")
  
  return(curationDrug)
}

curationDrug <- create_curationDrug(phenodata_DMData)

create_curationCell <- function(phenodata_DM, verbose = TRUE){
  curationCell <- unique(subset(phenodata_DM, select=c(cellid)))
  curationCell$dataset_cellid <- curationCell$cellid
  names(curationCell) <- c("unique.cellid", "dataset_cellid")
  rownames(curationCell) <- curationCell$unique.cellid
  
  return(curationCell)
}
curationCell <- create_curationCell(phenodata_DM)


create_curationTissue <- function(phenodata_DM, verbose = TRUE){
  curationTissue <- unique(subset(phenodata_DM, select=c(organ_id)))
  curationTissue$dataset_tissueid <- "Liver"
  names(curationTissue)[1] <- "unique.tissueid"
  rownames(curationTissue) <- "Hepatocyte"
  
  return(curationTissue)
}
curationTissue <- create_curationTissue(phenodata_DM)


#reading the metadata downloaded from diXa
s_Hepatocyte <- read.csv("../data/DM/s_Hepatocyte.csv", stringsAsFactors = F, sep = "\t")

create_drug <- function(metadataFile){
  
sub_hepa <- subset(metadataFile, subset = !duplicated(metadataFile$Factor.Value.Compound.),select = c("Factor.Value.Compound.","Term.Accession.Number.5"
                                                                                                      , "Characteristics.StdInChIKey.", "Comment.chEMBL.ID."), drop = F)
colnames(sub_hepa) <- c("dataset_drugid", "CHEBI.ID", "inchikey", "CHEMBL.ID")
#replace blank row with DMSO (important to add)
sub_hepa$dataset_drugid[which(sub_hepa$dataset_drugid == "")] <- "DMSO"
#merge to include unique drugid
mer_drug_cur <- merge(sub_hepa, curationDrug, by = "dataset_drugid", all.x = T)

colnames(mer_drug_cur)[5] <- "drugid" 
rownames(mer_drug_cur) <- mer_drug_cur$drugid
#reorder based on rows of curationDrug 
mer_drug_cur <- mer_drug_cur[rownames(curationDrug),]
return(mer_drug_cur)
}

drug <- create_drug(s_Hepatocyte)

create_cell <- function(phenodata_DM, verbose = TRUE){
  cell <- unique(subset(phenodata_DM,select=c(cellid, organ_id, species, test_type)))
  names(cell)<-c("cellid","tissueid", "species","testType")
  cell$tissueid<-"Liver"
  rownames(cell) <- cell$cellid
  
  return(cell)
}
cell <- create_cell(phenodata_DM)

drugMatrix <- ToxicoSet("drugMatrix_rat",
                  molecularProfiles=list("rna"= new_SE_DM),
                  cell=cell,
                  drug=drug,
                  sensitivityInfo=NULL,
                  sensitivityRaw=NULL,
                  sensitivityProfiles=NULL,
                  curationDrug=curationDrug,
                  curationCell=curationCell,
                  curationTissue=curationTissue,
                  datasetType = c("perturbation"),
                  verify = TRUE)

saveRDS(drugMatrix, "../results/drugMatrix.rds")



