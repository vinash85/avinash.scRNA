gs.enrichment.score = function(exp.count, genes.set){
	colSums(exp.count[rownames(exp.count) %in% genes.set, ])/colSums(exp.count)
}
norm.gs.enrichment.score = function(gs.score){
	(gs.score - min(gs.score))/(max(gs.score) - min(gs.score))
}



#' compute gamma delta T-cells score  
#' @param sco seurat object 
#' @examples
#' @export
get.gdT.score  = function(sco) {
	gdT.gs1 = c("CD3D", "CD3E", "TRDC", "TRGC1", "TRGC2")
	gdT.gs2= c("CD8A", "CD8B")
	sco$gdt.gs1  = gs.enrichment.score(sco@assays$RNA@counts, gdT.gs1) 

	sco$gdt.gs2 = (-1 * gs.enrichment.score(sco@assays$RNA@counts, gdT.gs2))

	sco$gdt = norm.gs.enrichment.score(sco$gdt.gs1) * norm.gs.enrichment.score(sco$gdt.gs2)
	sco
}

get.SingleR.annotation = function(sco) {
	require(SingleR)	
	singleR.annotations = list()
	hpca.se <- HumanPrimaryCellAtlasData()
	singleR.annotations$HumanPrimaryCellAtlasData.fine <- SingleR(test = sco@assays$RNA@data, ref = hpca.se, labels = hpca.se$label.fine)
	singleR.annotations$HumanPrimaryCellAtlasData.main <- SingleR(test = sco@assays$RNA@data, ref = hpca.se, labels = hpca.se$label.main)
	singleR.annotations
}

get.gdt.bloodatlas = function(sco){
	require(Matrix)
	gdt.blood.marker = fread("~/liulab_home/data/single_cell/markers/./blood_cell_category_rna_gdT-cell_Cell.tsv")
	exp.curr = sco@assays$RNA@data[rownames(sco) %in% gdt.blood.marker$Gene,]
tissue.specificity.matched = gdt.blood.marker[match(rownames(exp.curr),Gene)]$`RNA tissue specificity score` %>% as.numeric %>% ifelse(is.na(.), 1,.)
sco$gdt.specific.enhanced.exp = colSums(exp.curr, na.rm=T)
sco$gdt.specific.enhanced.exp.weighted = t(exp.curr) %*% tissue.specificity.matched %>% as.matrix %>% c
sco
}



get.TR.expression = function(sco, is.log=T){

	avg.expression <- function(sco, pattern, is.log=T){
		tryCatch({
			sel = grep(pattern, rownames(sco))
			aa = NULL

			if(length(sel) > 0){
				aa = sco@assays$RNA@data[sel,,drop=F]
				if(is.log) aa = exp(aa) - 1  
				if(ncol(aa) >1) aa= colSums(aa, na.rm=T)
			}
		unlist(aa)
	}, error = function(e) NULL) 
	} 
	for (gd in c("TRG", "TRD", "TRA", "TRB")){
		for (type in c("V", "D", "J", "C")){
			label = paste0(gd,type)
			out = avg.expression(sco, label, is.log=is.log)
			if(!is.null(out))
				sco %<>% AddMetaData(metadata=out,
					col.name =label)
		}
	}
	sco
}
