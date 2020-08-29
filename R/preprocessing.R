#' @export
preprocessing.scRNA <- function(counts, project.name="scRNA",  meta.data = NULL) {
    # require(Seurat)
    sco <- Seurat::CreateSeuratObject(counts = counts, project = project.name, min.cells = 3, min.features = 200,  meta.data = meta.data)
    # QC metrics 
    # We filter cells that have unique feature counts over 2,500 or less than 200
    # We filter cells that have >5% mitochondrial counts
    sco[["percent.mt"]] <- Seurat::PercentageFeatureSet(sco, pattern = "^MT-")
    # Seurat::VlnPlot(sco, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
    # Seurat::VlnPlot(sco, features = c( "percent.mt"))
    # sco <- subset(sco, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
    # normalizing dataset 
    sco <- Seurat::NormalizeData(sco, normalization.method = "LogNormalize", scale.factor = 10000)
    all.genes <- rownames(sco)
    sco <- Seurat::ScaleData(sco, features = all.genes)
    sco 
}

## Linear mixed model with NMLE faster than LME4

#' @export
extract_nmle_table <- function (m1){
    mod = summary(m1)
    beta <- m1$coefficients$fixed[2] #$fixed is not needed
    se <- m1$varFix[2]
    t <- beta/se
    p<- anova(m1)[2,4]
    table=data.frame(cbind(beta,se,t,p))
    return(table)
}
#' @export
eval.nlme = function(gene.expression, data.dt){
library(nlme)
    # data.dt = data.table with colums  Response1 and dataset (factor to control)
    gene.expression = avinash::qnorm.array(gene.expression)
    tryCatch(
        {
            data.dt$col = gene.expression
            data.dt = data.dt[!is.na(col)]
            data.dt = data.dt[dataset %in% names(which(table(data.dt$dataset) > 10))]  # only choose those with at least 10 obervations 
            m1 <- nlme::lme(col~ Response1, random=~1|dataset, data=data.dt)
            extract_nmle_table(m1)
        },error=  function(e) rep(NA,4))
}


preprocessing.standard.seurat <- function(sco, num.dim=50) {

  sco %<>% Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(verbose = FALSE) 
  sco %>% RunPCA(pc.genes = sco@meta.data, npcs = num.dim, verbose = FALSE)%>%
  RunUMAP(reduction = "pca", dims = seq(num.dim)) %>% 
  RunTSNE(reduction = "pca", dims = seq(num.dim)) %>% 
  FindNeighbors(reduction = "pca", dims = seq(num.dim)) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()
}


preprocessing.harmony.seurat <- function(sco, num.dim=50, patient.name="patient.name"){
    require(Seurat)
    require(magrittr)
    require(harmony)
    sco[["percent.mt"]] <- PercentageFeatureSet(sco, pattern = "^MT-")
    sco %>%   subset(subset = percent.mt >= 4) %>%
    Seurat::NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
    ScaleData(verbose = FALSE) %>% 
    RunPCA(pc.genes = sco@var.genes, npcs = num.dim, verbose = FALSE)%>%
    RunHarmony(patient.name) %>%
    RunUMAP(reduction = "harmony", dims = seq(num.dim)) %>% 
    RunTSNE(reduction = "harmony", dims = seq(num.dim)) %>% 
    FindNeighbors(reduction = "harmony", dims = seq(num.dim)) %>% 
    FindClusters(resolution = 0.5) %>% 
    identity()
}