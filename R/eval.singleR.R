myFastCor.multicores = function(x,y, method="spearman", num=20000, nthreads=1){
    if(is.null(num))  num=0
    if(nthreads>1){
      split.size = floor(ncol(y)/nthreads)
        splits = lapply(seq(nthreads), function(tt){
          start = (tt-1)* split.size +1
          end = ifelse(tt< nthreads, tt*split.size, ncol(y))
          seq(start, end)
        })
        out = mclapply(splits, function(tt) myFastCor(x=x,y=y[,tt], method=method, num=num), mc.cores = nthreads) %>%
            do.call(cbind,.) %>% set_colnames(colnames(y))
    }else{
        out = myFastCor(x=x,y=y, method=method, num=num)
    }
    out
}


myFastCor <- function(x,y, method="spearman", num=20000) {
    num = max(c(num, 10))
    if (ncol(x) < num | ncol(y) < num){
        out = cor(x, y, method="spearman", use="pairwise.complete.obs")
    }else{
        if(method=="spearman"){
            x%<>%
                apply(.,2, rank, na.last="keep")
            y%<>%
                apply(.,2, rank, na.last="keep")
        }
        out = WGCNA::cor(x, y, use="pairwise.complete.obs")
    }
    out
}
calcCorEnrich <- function(x,y, num=NULL, nthreads=1) {
  rownames(x) %<>% toupper
  rownames(y) %<>% toupper
    common = intersect(rownames(x), rownames(y))
    myFastCor.multicores(x[common,, drop=F], y[common,,drop=F], num=num, nthreads = nthreads)
}



mytheme <- function() {
    theme(panel.border = element_blank(), panel.background = element_blank(), panel.grid.major = element_line(colour = "grey"),
      panel.grid.minor = element_line(colour = "grey"), axis.line = element_line(colour = "black"))
}
combined.singleR.cor.datasets =  function(singleR.list) {
    enrich.dt = do.call(rbind, singleR.list)
    enrich.dt[,label.main:=gsub("_cell", label.main, replacement = "-cell", ignore.case = T)]
    enrich.dt[,label.main:=gsub(" cell", label.main, replacement = "-cell", ignore.case = T)]
    enrich.dt[,label.main:=gsub("cells", label.main, replacement = "cell", ignore.case = T)]
    enrich.dt[,label.fine:=gsub("_cell", label.fine, replacement = "-cell", ignore.case = T)]
    enrich.dt[,label.fine:=gsub(" cell", label.fine, replacement = "-cell", ignore.case = T)]
    enrich.dt[,label.fine:=gsub("cells", label.fine, replacement = "cell", ignore.case = T)]
    enrich.dt
}

plot.mat.singleR.enrichment <- function(mat, data.name, dirname, prefix="", width=20, height = 16, dt2  = NULL) {
    # ref.se <- BlueprintEncodeData(rm.NA = "none")
    require(SingleR)
    require(ggplot2)
    dir.create(dirname, recursive = T)
    prefix = sprintf("%s/%s", dirname, prefix)
    if(is.null(dt2)){
        ref.se =  eval(parse(text=paste0(data.name, "()"))) 
        col.dat = colData(ref.se)
        col.dat$Var1= paste(rownames(col.dat),seq(nrow(col.dat)))  
        col.dat = as.data.table(col.dat)
        aa = assay(ref.se) %>%
        t() %>% scale(., center=T, scale=T) %>% t()
        rownames(aa) %<>% toupper
        colnames(aa) = col.dat$Var1
        
        # aa = aa[rownames(aa)%in%toupper(VariableFeatures(lee.sco)),]
        
        mat.signature.enrich = calcCorEnrich(aa, mat)
        temp = melt(mat.signature.enrich) 
        dt2= col.dat[match(as.character(temp$Var1),Var1), .(label.main,label.fine)] %>%
        cbind(., temp) 
    }
    # mat.enrich = mat.signature.enrich[,"combined"]
    # dt1 = cbind(col.dat, mat.enrich) %>% as.data.table()
    
    dt.sum = dt2[,.(enrich.mean=mean(value)), by=label.main] %>%
    .[order(enrich.mean, decreasing =T)]
    dt.sum2 = dt2[,.(enrich.mean=mean(value)), by=label.fine] %>%
    .[order(enrich.mean, decreasing =T)]
    dt2$label.main = factor(dt2$label.main, level=dt.sum$label.main)
    dt2$label.fine = factor(dt2$label.fine, level=dt.sum2$label.fine)
    p3= ggplot(dt2, aes(y=value,x=label.main))  +   
    geom_boxplot(aes(fill=label.main)) + facet_wrap(~Var2, scales = "free") + 
    theme(axis.text.x = element_blank()) +
    mytheme() 
    ggsave(p3,filename = sprintf("%s_cor3.main.pdf", prefix),  width=width, height = height)
    
    p3= ggplot(dt2, aes(y=value,x=label.main))  +   
    geom_boxplot(aes(fill=Var2)) + 
    theme(axis.text.x = element_blank()) +
    mytheme() + 
    theme(axis.text.x = element_text(angle=60, size = 8))
    ggsave(plot = p3,filename = sprintf("%s_cor2.main.pdf", prefix),  width=width, height = height)
    
    p4= ggplot(dt2, aes(y=value,x=label.fine))  +   
    geom_boxplot(aes(fill=label.main)) + 
    theme(axis.text.x=element_blank()) +
    mytheme() + 
    facet_wrap(~Var2) + theme(axis.text.x = element_text(angle=60, size = 5)) 
    ggsave(p4,filename = sprintf("%s_cor3.fine.pdf", prefix),  width=2*width, height = 2*height, limitsize = F)
    
    p4= ggplot(dt2, aes(y=value,x=label.fine))  +   
    geom_boxplot(aes(fill=Var2)) + 
    theme(axis.text.x=element_blank()) +
    mytheme() + 
    theme(axis.text.x = element_text(angle=60, size = 5))
    ggsave(p4,filename = sprintf("%s_cor2.fine.pdf", prefix),  width=width, height = height, limitsize = F)

    dt2
}

#'  SingleR analysis from Differential expression 
#' @param dataset.deg.mat a matrix of differential expression with genes (rows) columns (contrasts)
#' @param dirname dirname to stor the result 
#' @return a data table with combined analysis  
#' @export
#' @examples
eval.singleR.deg = function(dataset.deg.mat, dirname, width=16, height=10){
    quantile.signed = function(val, probs){
        m = mean(val)
        if(m<0) probs = 1-probs
        quantile(val, probs)
    }

    library(SingleR)
    library(ggrepel)
    rownames(dataset.deg.mat) %<>% toupper
    singleR.data.names = grep("Data$", ls("package:SingleR"), value=T)
    dataset.singleR.list = list()
    for (ii in singleR.data.names) {
        dataset.singleR.list[[ii]] = plot.mat.singleR.enrichment(mat=dataset.deg.mat, data.name=ii, dirname= dirname, prefix = ii)
    }


    dataset.enrich.dt = combined.singleR.cor.datasets(dataset.singleR.list)

    aa = plot.mat.singleR.enrichment(mat=NULL, data.name="temp", dirname=dirname, prefix = "singleRcombined", width=width, height= height, dt2 = dataset.enrich.dt)

    combined.dataset.enrich.sum = dataset.enrich.dt[, .(enrich.mean=mean(value),
        enrich.median=quantile(value,probs=0.5),
        enrich.75=quantile(value,probs=0.75),
        enrich.75.signed=quantile.signed(value,probs=0.75),
        enrich.25=quantile(value,probs=0.25),
        label.main=unique(label.main)
        ), by=label.fine]

    dataset.merged = combined.dataset.enrich.sum
    dataset.merged.sub = dataset.merged[enrich.mean > 0.1 | enrich.mean < -0.075| abs(enrich.median) > 0.1]
    p1 = ggplot(dataset.merged, aes(y=enrich.mean,x=enrich.median))  +   
    geom_point() + 
    geom_text_repel(
        data = dataset.merged.sub,
        aes(label =label.fine, color = factor(label.main)),
        size = 3,
        box.padding = unit(0.35, "lines"),
        point.padding = unit(0.3, "lines")
        ) +
    theme_bw() 
    ggsave(sprintf("%s/dataset.merged.pdf", dirname), p1, width=width, height = height)
    dataset.merged

}
