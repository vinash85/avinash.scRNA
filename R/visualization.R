
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




