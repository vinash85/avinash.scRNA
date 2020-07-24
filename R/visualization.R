
#' compute t-SNE or umap of each round of HIPPO
#' @param sce SingleCellExperiment object with hippo object in it.
#' @param method a string that determines the method for dimension
#' reduction: either 'umap' or 'tsne
#' @param perplexity numeric perplexity parameter for Rtsne function
#' @param featurelevel the round of clustering that you will extract
#' features to reduce the dimension
#' @return a data frame of dimension reduction result for each
#' k in 1, ..., K
#' @examples
#' data(toydata)
#' set.seed(20200321)
#' set.seed(20200321)
#' toydata = hippo(toydata,K = 10,z_threshold = 1,outlier_proportion = 0.01)
#' toydata = hippo_dimension_reduction(toydata, method="tsne")
#' hippo_tsne_plot(toydata)
#' @export
myhippo_dimension_reduction = function(sce, method = c("umap", "tsne"),
 perplexity = 30,
 featurelevel = 1, n_threads = 32,
 ...) {
  hippo_object = sce@int_metadata$hippo
  dflist = list()
  K = ncol(hippo_object$labelmatrix)
  if (method == "umap"){
    dimred = log(as.matrix( hippo_object$X[hippo_object$features[[1]]$gene,] ) + 1) 
    dimred = uwot::umap(t(dimred),n_threads=n_threads, ...)
    # uwot::umap(.,n_threads=32, ...)
  }else{
    dimred = tsne = Rtsne::Rtsne(log(t(hippo_object$X[hippo_object$features[[1]]$gene,
      ]) + 1), perplexity = perplexity,
    check_duplicates = FALSE)$Y
  }
  dimred = as.data.frame(dimred)
  dimreddf = data.frame()
  for (i in 2:K){
    df = preprocess_homogeneous(sce, label = hippo_object$labelmatrix[,i])
    df$selected_feature = df$gene %in% hippo_object$features[[i -1]]
    df$K = i
    dflist[[i]] = df
    dimreddf = rbind(dimreddf,
     data.frame(dim1 = dimred$V1, dim2 = dimred$V2,
      K = i,
      label = hippo_object$labelmatrix[, i]))
  }
  if (method == "umap"){
    sce@int_metadata$hippo$umap = NA
    colnames(dimreddf) = c("umap1", "umap2", "K", "label")
    dimreddf$label = as.factor(dimreddf$label)
    sce@int_metadata$hippo$umap = dimreddf
  }else{
    sce@int_metadata$hippo$tsne = NA
    colnames(dimreddf) = c("tsne1", "tsne2", "K", "label")
    dimreddf$label = as.factor(dimreddf$label)
    sce@int_metadata$hippo$tsne = dimreddf
  }
  return(sce)
}


myone_level_clustering = function(subX, z_threshold) {
  subdf = preprocess_heterogeneous(subX)
  subdf = compute_test_statistic(subdf)
  features = subdf[subdf$zvalue > z_threshold, ]
  nullfeatures = data.frame(matrix(ncol = 11, nrow = 0))
  colnames(nullfeatures) = c("gene", "gene_mean", "zero_proportion",
   "gene_var", "samplesize", "expected_pi", "se",
   "minus_logp","zvalue", "subsetK", "K")
  if (nrow(features) < 10) {
    return(list(features = nullfeatures, pcs = NA, km = NA))
  }
  if (nrow(features) < 10) {
    return(list(features = nullfeatures, pcs = NA, km = NA,
      unscaled_pcs = NA,subdf = NA))
  }
  pcs = tryCatch(expr = {
    irlba::irlba(log(subX[features$gene, ] + 1), min(9, nrow(features) -
     1, ncol(subX) - 1))$v
  }, error = function(e) {
    NA
  }, warning = function(w) {
    NA
  })
  if (is.na(pcs[1])) {
    return(list(features = nullfeatures, pcs = NA, km = NA,
      unscaled_pcs = NA,
      subdf = NA))
  } else {
    unscaledpc = irlba::prcomp_irlba(log(Matrix::t((subX[features$gene,])) + 1),
     n = min(9, nrow(features) - 1,
       ncol(subX) - 1),
     scale. = FALSE, center = FALSE)$x
    km = kmeans(pcs, 2, nstart = 10, iter.max = 50)
  }
  return(list(features = features,
    pcs = pcs,
    km = km,
    unscaled_pcs = unscaledpc,
    subdf = subdf))
}



#' HIPPO's hierarchical clustering
#'
#' @param sce SingleCellExperiment object
#' @param K number of clusters to ultimately get
#' @param z_threshold numeric > 0 as a z-value threshold
#' for selecting the features
#' @param outlier_proportion numeric between 0 and 1, a cut-off
#' so that when the proportion of important features reach this
#' number, the clustering terminates
#' @param verbose if set to TRUE, it shows progress of the algorithm
#' @examples
#' data(toydata)
#' toydata = hippo(toydata,K = 10,z_threshold = 1,outlier_proportion = 0.01)
#' @return a list of clustering result for each level of k=1, 2, ... K.
#' @export
myhippo = function(sce, K = 20,
 z_threshold = 2,
 outlier_proportion = 0.001,
 verbose = TRUE) {
  if (is(sce, "SingleCellExperiment")) {
    X = sce@assays@data$counts
  } else if (is(sce, "mat
    rix")) {
    sce = SingleCellExperiment::SingleCellExperiment(assays=list(counts = sce))
    X = sce@assays@data$counts
  } else {
    stop("input must be either matrix or SingleCellExperiment object")
  }
  if (outlier_proportion > 1 | outlier_proportion < 0) {
    stop("Outlier_proportion must be a number between 0 and 1.
     Default is 5%")
  }
  param = list(z_threshold = z_threshold,
   outlier_proportion = outlier_proportion,
   maxK = K)
  outlier_number = nrow(X) * outlier_proportion
  labelmatrix = matrix(NA, ncol(X), K)
  labelmatrix[, 1] = 1
  eachlevel = list()
  subX = X
  subXind = seq(ncol(X))
  withinss = rep(0, K)
  oldk = 1
  features = list()
  featuredata = list()
  for (k in 2:K) {
    thisk = myone_level_clustering(subX, z_threshold)
    if (is.na(thisk$features$gene[1])) {
      if(verbose){
        message("not enough important features left; terminate the procedure")
      }
      labelmatrix = labelmatrix[, seq((k - 1))]
      break
    }
    if (nrow(thisk$features) < outlier_number) {
      if(verbose){
        message("not enough important features; terminate the procedure")
      }
      labelmatrix = labelmatrix[, seq((k - 1))]
      break
    }
    # if(min(table(thisk$km$cluster)) <= 1){
    #   if(verbose){
    #     message("only one cell in a cluster: terminate procedure")
    #   }
    #   labelmatrix = labelmatrix[, seq((k - 1))]
    #   break
    # }
    if (verbose) {message(paste0("K = ", k, ".."))}
    labelmatrix[, k] = labelmatrix[, k - 1]
    labelmatrix[subXind[thisk$km$cluster == 2], k] = k
    oneind = thisk$km$cluster == 1
    twoind = thisk$km$cluster == 2
    if(sum(oneind) >= 2){
      withinss[oldk] = sum(apply(thisk$unscaled_pcs[oneind, ],1, var)^2)
    }else{
      withinss[oldk] = var(as.numeric(thisk$unscaled_pcs[oneind,]))^2
    }
    if(sum(twoind) >= 2){
      withinss[k] = sum(apply(thisk$unscaled_pcs[twoind, ], 1, var)^2)
    }else{
      withinss[k] = var(as.numeric(thisk$unscaled_pcs[twoind,]))^2
    }


    ind = which(table(thisk$km$cluster) <= 5)
    if (length(ind) >= 1){
      valid_indices = seq(k-1)[-ind]
      oldk = which(withinss == max(withinss[valid_indices]))
    }else{
      oldk = which.max(withinss[seq(k-1)])
    }
    if (sum(labelmatrix[, k] == oldk) < 2) {
      if(verbose){
        message("too few cells in one cluster; terminating the procedure")
      }
      labelmatrix = labelmatrix[, seq(k)]
      break
    }
    subX = X[thisk$features$gene, which(labelmatrix[, k] == oldk)]
    subXind = which(labelmatrix[, k] == oldk)
    thisk$features$subsetK = oldk
    thisk$features$K = k
    features[[k - 1]] = thisk$features
  }
  sce@int_metadata$hippo = list(X = X,
    features = features,labelmatrix = labelmatrix,
    z_threshold = z_threshold, param = param,
    outlier_proportion = outlier_proportion)
  return(sce)
}



#' compute tile contour tanaka plots from cordinates 
#' @param x  x cooridnates 
#' @param y  y cooridnates 
#' @param z  value to show 
#' 
#' @return a list with 
#' predicted data.table 
#' tile, conotour and tananka plot 
#' @export
plot.contour = function(x, y, z, cutoff = 0.1, fmla= as.formula(z ~1), interpolator="gam", grid.size=100, nbins=10, bw.scale=1){

  library(ggplot2)
  library(scales)
  library(directlabels)
  library(metR)
  dat = data.frame(x=x,y=y,z=z)
  bw.scale.fn = function(x)  MASS::bandwidth.nrd(x) * bw.scale
  h.scale = apply(dat,2,bw.scale.fn)
  kd.out = MASS::kde2d(dat$x, dat$y, n=c(grid.size,grid.size), h=h.scale)

  ran = range(c(kd.out$z))
  cutoff.relax = ran[1] + (ran[2] - ran[1]) * cutoff/10 
  cutoff = ran[1] + (ran[2] - ran[1]) * cutoff 
  # seq(cutoff[1],cutoff[2], length.out=20)[3]
  newdata = expand.grid(kd.out$x, kd.out$y) %>% 
  set_colnames(c("x", "y")) %>% as.data.table %>%
  .[,density:=c(kd.out$z)] 

  if(interpolator=="gam"){
    predicted = interp.gam(dat, newdata)
  }else if(interpolator=="kriged"){
    predicted = interp.kriged(dat, newdata)
  }else if(interpolator=="kriged.variogram"){
    predicted = interp.kriged.variogram(dat, newdata, fmla=fmla)
  }
  # browser()
  out = newdata %>% as.data.table %>%
  .[,var1.pred:=predicted] %>% 
  .[,var1.pred.round:=round(var1.pred, digit=3)] %>%
  .[,var1.norm:=scale(var1.pred, scale=T, center=T)]

  p1 = out[density > cutoff] %>%
  ggplot(aes(x=x, y=y)) + geom_tile(aes(fill=var1.pred)) +
   # coord_equal() +
  scale_fill_gradient(low = "yellow", high="red") +
  scale_x_continuous(labels=comma) + scale_y_continuous(labels=comma) +
  theme_bw()

  breaks = round(quantile(out$var1.pred, seq(0, 1, length.out=nbins)), 2)

  p2 = out %>%
  ggplot( aes(x=x, y=y)) + 
  stat_contour(aes(z=var1.pred.round, colour = ..level..), bins = nbins) + 
  scale_fill_gradient(low = "yellow", high="red") +
  scale_x_continuous(labels=comma) + scale_y_continuous(labels=comma) +
  stat_subset(aes(subset = density < cutoff), geom = "raster", 
    fill = "#EBEBEB", alpha=0.85) +
   # geom_text_contour(aes(z=var1.pred.round, color=NULL), bins = nbins, skip=2, size=1) +
  theme_bw()
  # p3 =  direct.label(p2, "bottom.pieces")

  
  p4= ggplot( out, aes(x, y, z = var1.norm)) +
  geom_contour_fill(na.fill = TRUE, bins=nbins) +
  geom_contour_tanaka(bins=nbins) +
  scale_fill_divergent() +
  # scale_x_longitude() +
  # scale_y_latitude()   + 
  stat_subset(aes(subset = density < cutoff), geom = "raster", 
    fill = "#EBEBEB") 


  list(predicted.dt = out, tile.plot = p1, contour.plot=p2,tanaka.plot=p4, cutoff=cutoff, cutoff.relax=cutoff.relax)
}


#' interploate using kriging variagram  
#' @param x  x cooridnates 
#' @param y  y cooridnates 
#' @param z  value to show 
#' 
#' @return predicted data.table 
interp.kriged.variogram = function(dat, newdata, fmla = as.formula(z ~1)){
  library(metR)
  library(sp)
  library(gstat)
  library(magrittr)
  coordinates(dat) <- ~ x + y
  dat %<>% cbind(., dat@coords)
  coordinates(newdata) <- ~ x + y
  newdata %<>% cbind(., newdata@coords)
  lzn.vgm <- variogram(fmla, dat) # calculates sample variogram values 
  # lzn.fit <- fit.variogram(lzn.vgm, model=vgm(1, "Sph", 900, 1)) # fit model
  options(warn = 0)
  lzn.fit <- fit.variogram(lzn.vgm, model=vgm(c("Exp", "Ste", "Sph")),  fit.kappa = TRUE) # fit model
  if(nrow(dat) > 3000){
    dat = dat[sample(nrow(dat), 3000),]
  }

  lzn.kriged <- krige(fmla, dat, newdata, model=lzn.fit)
  lzn.kriged$var1.pred 
}

interp.kriged = function(dat, newdata){
  epsilon <- 1e-3
  require(DiceKriging)
# fit a kriging model 
# m.1 <- km(~.^2,design=dat[,c("x","y")], response=dat$z)
  m.1 <- km(design=dat[,c("x","y")], response=dat$z)
# estimate a response 
  aa = predict(m.1, newdata=newdata[,c("x","y")], type="UK",    se.compute=FALSE, nugget=epsilon)
  aa$mean

}



interp.gam = function(dat, newdata){
  require(mgcv)
# + te(x,y,bs = rep("cr",2),fx=TRUE)
  foo <- gam(z~s(x,bs="cc",fx=TRUE)+
    s(y,bs="cc",fx=TRUE), 
  # te(x,y,bs=rep("cr",2)),
  # s(x*y,bs="cr",fx=TRUE)+
  # s(x*y,bs="cr",fx=TRUE),
    data=dat) 
  predict(foo,newdata)
}


ContourDimPlot = function(
  object,
  features,
  object.contour = NULL,
  dims = c(1, 2),
  cells = NULL,
  cols = NULL,
  pt.size = NULL,
  reduction = NULL,
  group.by = NULL,
  split.by = NULL,
  shape.by = NULL,
  interpolator = "gam",
  order = NULL,
  label = FALSE,
  label.size = 4,
  slot = 'data',
  repel = FALSE,
  cells.highlight = NULL,
  cols.highlight = '#DE2D26',
  sizes.highlight = 1,
  na.value = 'grey50',
  combine = TRUE,
  ncol = NULL,
  ...
  ) {
  require(purrr)
  if (length(x = dims) != 2) {
    stop("'dims' must be a two-length vector")
  }
  p =DimPlot(object = object, reduction = reduction, label = TRUE, pt.size = 0.5) + NoLegend()
  if(!is.null(object.contour)) object = object.contour
  reduction <- reduction %||% DefaultDimReduc(object = object)
  cells <- cells %||% colnames(x = object)
  data <- Embeddings(object = object[[reduction]])[cells, dims]
  data <- as.data.frame(x = data)
  dims <- paste0(Key(object = object[[reduction]]), dims)

  data <- FetchData(
    object = object,
    vars = c(dims, features),
    cells = cells,
    slot = slot
    )

  contour.out = plot.contour(data[[dims[[1]]]], data[[dims[[2]]]], data[[features]], interpolator=interpolator,  cutoff=0.3, grid.size=100, nbins=20, bw.scale=0.6)
  
  layers = lapply(p$layers, function(tt) {
    curr.class =  tt$data %>% class
    if(curr.class=="waiver")
      tt$data = p$data
    tt
  })

  contour.out.curr = contour.out$predicted.dt
  layers[[1]]$aes_params$alpha = .4
  p1 =  ggplot()+
  geom_contour_fill(data = contour.out.curr, aes(x, y, z=var1.norm), na.fill = TRUE, bins=20) +
  geom_contour_tanaka(data = contour.out.curr, aes(x, y, z=var1.norm), bins=20) +
  scale_fill_divergent() + layers + theme_bw()
  return(p1) 

} 