#' Differential expression using DESeq2
#'
#' Identifies differentially expressed genes between two groups of cells using
#' DESeq2
#'
#' @references Love MI, Huber W and Anders S (2014). "Moderated estimation of
#' fold change and dispersion for RNA-seq data with DESeq2." Genome Biology.
#' https://bioconductor.org/packages/release/bioc/html/DESeq2.html
#' @param data.use Data matrix to test
#' @param cells.1 Group 1 cells
#' @param cells.2 Group 2 cells
#' @param verbose Print a progress bar
#' @param ... Extra parameters to pass to DESeq2::results
#' @return Returns a p-value ranked matrix of putative differentially expressed
#' genes.
#'
#' @details
#' This test does not support pre-filtering of genes based on average difference
#' (or percent detection rate) between cell groups. However, genes may be
#' pre-filtered based on their minimum detection rate (min.pct) across both cell
#' groups. To use this method, please install DESeq2, using the instructions at
#'  https://bioconductor.org/packages/release/bioc/html/DESeq2.html
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   pbmc_small
#'   DESeq2DETest(pbmc_small, cells.1 = WhichCells(object = pbmc_small, idents = 1),
#'               cells.2 = WhichCells(object = pbmc_small, idents = 2))
#' }

DESeq2DETest <- function(
  data.use,
  cells.1,
  cells.2,
  verbose = TRUE,
  ...
) {
  # if (!PackageCheck('DESeq2', error = FALSE)) {
  #   stop("Please install DESeq2 - learn more at https://bioconductor.org/packages/release/bioc/html/DESeq2.html")
  # }
  group.info <- data.frame(row.names = c(cells.1, cells.2))
  group.info[cells.1, "group"] <- "Group1"
  group.info[cells.2, "group"] <- "Group2"
  group.info[, "group"] <- factor(x = group.info[, "group"])
  group.info$wellKey <- rownames(x = group.info)
  dds1 <- DESeq2::DESeqDataSetFromMatrix(
    countData = data.use,
    colData = group.info,
    design = ~ group
  )
  dds1 <- DESeq2::estimateSizeFactors(object = dds1)
  dds1 <- DESeq2::estimateDispersions(object = dds1, fitType = "local")
  dds1 <- DESeq2::nbinomWaldTest(object = dds1)
  res <- DESeq2::results(
    object = dds1,
    contrast = c("group", "Group1", "Group2"),
    alpha = 0.05,
    ...
  )
  # to.return <- data.frame(p_val = res$pvalue, row.names = rownames(res))
  return(res)
}

#'  Differential expression using DESeq2 using Seurat package
#' @param sco Seurat object. With response variable in "CR", "PR", "PD" and "SD"
#' @param method a string that determines the method for dimension
#' @return a data table with Seurat output 
#' @export
#' @examples
#' toydata = mydeg(sco) 
mydeg <- function(sco, 
  responder.string=c("CR", "PR","responder", "responders"), 
  nonresponder.string=c("PD", "SD", "Non-Responder", "Non-Responders") 
  ) {
  responder.string = toupper(responder.string)
  nonresponder.string = toupper(nonresponder.string)
  temp= toupper(sco$response) %>%
  gsub(" ", ., replacement="")
  sco$response = temp
  valid.string = c(responder.string, nonresponder.string)
  sel.inx = which(sco$response %in% c(responder.string, nonresponder.string))
  exp.curr1 = sco@assays$RNA@counts
  meta.dt1 = sco@meta.data %>%
  as.data.table()
  if(length(sel.inx) < ncol(sco)){
    exp.curr1 = exp.curr1[, sel.inx]
    meta.dt1=meta.dt1[sel.inx,]
  }
  meta.dt1=meta.dt1[,.(binaryResponse=ifelse(response %in% responder.string,1 ,0) , patient=patient.name)] 

  meta.curr = list()
  exp.curr2 = list()
  for(patient in unique(meta.dt1$patient)){
    inx = which(meta.dt1$patient==patient)
    if(length(inx)> 1){
      exp.curr2[[patient]] = rowSums(exp.curr1[,inx],na.rm=T)
    }else{
      exp.curr2[[patient]] = exp.curr1[,inx]
    }
    meta.curr[[patient]] = meta.dt1[inx[1],]
  }
  meta.dt = do.call(rbind, meta.curr)
  exp.curr = do.call(cbind, exp.curr2) %>% round
  responders = meta.dt[binaryResponse==1]$patient %>% as.character
  nonresponders = meta.dt[binaryResponse==0]$patient%>% as.character
  deseq.out = DESeq2DETest(data.use=exp.curr[,c(responders,nonresponders)], cells.1=responders, cells.2=nonresponders)
  deseq.dt = deseq.out %>%
  as.data.frame() %>%
  mutate(gene=rownames(.)) %>%
  data.table() %>% 
  .[order(pvalue)]
  deseq.dt
}

