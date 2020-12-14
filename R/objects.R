# Call gc() to perform garbage collection
#
CheckGC <- function() {
  if (getOption(x = "Seurat.memsafe")) {
    gc(verbose = FALSE)
  }
}

#' @importFrom methods is
#' @importFrom Matrix sparseMatrix
#' @importFrom utils packageVersion
#'
#' @rdname h5ad
#' @export
#' @method ReadH5AD H5File
#'
ReadH5AD.new <- function(file, assay = 'RNA', verbose = TRUE, ...) {
  # Pull assay data
  # If X is an H5D, assume scaled
  # Otherwise, if file$exists(name = 'raw'), assume X is normalized
  # Otherwise, assume file[['X']] is raw counts
  file <- hdf5r::h5file(filename = file, mode = 'r')
  if (verbose) {
    message("Pulling expression matrices and metadata")
  }
  if (is(object = file[['X']], class2 = 'H5Group')) {
    x <- as.sparse(x = file[['X']])
  } else {
    x <- file[['X']][, ]
  }
  # x will be an S3 matrix if X was scaled, otherwise will be a dgCMatrix
  scaled <- is.matrix(x = x)
  if (verbose) {
    message("Data is ", ifelse(test = scaled, yes = 'scaled', no = 'unscaled'))
  }
  # Pull cell- and feature-level metadata
  obs <- file[['obs']][]
  x.var <- file[['var']][]
  x.var$index.new =  x.var$index
  inx = duplicated(x.var$index.new)
  x.var$index.new[inx] %<>% paste(.,"_1") 
  rownames(x = x) <- rownames(x = x.var) <- x.var$index.new
  colnames(x = x) <- rownames(x = obs) <- obs$index
  # Pull raw expression matrix and feature-level metadata
  if (file$exists(name = 'raw.X')) {
    if (is(object = file[['raw.X']], class2 = 'H5Group')) {
      raw <- as.sparse(x = file[['raw.X']])
    } else {
      raw <- as.sparse(file[['raw.X']][, ])
    }

    #raw <- as.sparse(x = file[['raw.X']])
    raw.var <- file[['raw.var']][]
    slot(object = raw, name = 'Dim') <- c(nrow(x = raw.var), nrow(x = obs))
      raw.var$index.new =  raw.var$index
  inx = duplicated(raw.var$index.new)
    raw.var$index.new[inx] %<>% paste(.,"_1") 
    rownames(x = raw) <- rownames(x = raw.var) <- raw.var$index.new
    colnames(x = raw) <- obs$index
    raw.var <- raw.var[, -which(x = colnames(x = raw.var) == 'index'), drop = FALSE]
    x.slot <- ifelse(test = scaled, yes = 'scale.data', no = 'data')
  } else {
    # If X is scaled, we required normalized data present in raw
    if (scaled) {
      stop("Seurat requires normalized data present in the raw slot when X is scaled")
    } else {
      x.slot <- 'raw'
    }
  }
  obs <- obs[, -which(x = colnames(x = obs) == 'index'), drop = FALSE]
  x.var <- x.var[, -which(x = colnames(x = x.var) == 'index'), drop = FALSE]
  # Merge raw.var and x.var
  # Only happens when we have a raw.X and raw.var in the h5ad file
  if (x.slot != 'raw') {
    if (verbose) {
      message("Merging feature-level metadata dataframes")
    }
    x.var <- x.var[, -which(x = colnames(x = x.var) %in% colnames(x = raw.var))]
    meta.features <- merge(x = raw.var, y = x.var, by = 0, all = TRUE)
    rownames(x = meta.features) <- meta.features$Row.names
    meta.features <- meta.features[, -which(x = colnames(x = meta.features) == 'Row.names'), drop = FALSE]
    rm(raw.var)
  } else {
    meta.features <- x.var
  }
  # Fix meta feature colnames
  colnames(x = meta.features) <- gsub(
    pattern = 'dispersions_norm',
    replacement = 'dispersion.scaled',
    x = colnames(x = meta.features)
  )
  colnames(x = meta.features) <- gsub(
    pattern = 'dispersions',
    replacement = 'dispersion',
    x = colnames(x = meta.features)
  )
  colnames(x = meta.features) <- gsub(
    pattern = 'means',
    replacement = 'mean',
    x = colnames(x = meta.features)
  )
  colnames(x = meta.features) <- gsub(
    pattern = '_',
    replacement = '.',
    x = colnames(x = meta.features)
  )
  if ('highly.variable' %in% colnames(x = meta.features)) {
    meta.features$highly.variable[is.na(x = meta.features$highly.variable)] <- FALSE
  }
  rm(x.var)
  CheckGC()
  # Fix metadata colnames
  colnames(x = obs) <- gsub(
    pattern = '_',
    replacement = '.',
    x = colnames(x = obs)
  )
  colnames(x = obs) <- gsub(
    pattern = 'n.genes',
    replacement = paste0('nFeatures_', assay),
    x = colnames(x = obs)
  )
  colnames(x = obs) <- gsub(
    pattern = 'n.counts',
    replacement = paste0('nCount_', assay),
    x = colnames(x = obs)
  )
  # Assemble assay object
  if (verbose) {
    message("Creating assay object")
    message(
      "Storing X as ",
      x.slot,
      ifelse(
        test = x.slot != 'counts',
        yes = paste(" and raw as", ifelse(test = scaled, yes = 'data', no = 'counts')),
        no = ''
      )
    )
  }
  if (scaled) {
    assays <- list(CreateAssayObject(data = raw))
    assays[[1]] <- SetAssayData(
      object = assays[[1]],
      slot = 'scale.data',
      new.data = x
    )
    rm(raw)
  } else if (x.slot == 'data') {
    assays <- list(CreateAssayObject(counts = raw))
    assays[[1]] <- SetAssayData(
      object = assays[[1]],
      slot = 'data',
      new.data = x
    )
    rm(raw)
  } else {
    assays <- list(CreateAssayObject(counts = x))
  }
  names(x = assays) <- assay
  # Add meta feature information
  #assays[[assay]][[names(x = meta.features)]] <- meta.features
  # Add highly variable feature information
  if ('highly.variable' %in% colnames(x = assays[[assay]][[]])) {
    if (verbose) {
      message("Setting highly variable features")
    }
    hvf.info <- HVFInfo(object = assays[[assay]])
    hvf.info <- hvf.info[order(hvf.info$dispersion, decreasing = TRUE), , drop = FALSE]
    means.use <- (hvf.info$mean > 0.1) & (hvf.info$mean < 8)
    dispersions.use <- (hvf.info$dispersion.scaled > 1) & (hvf.info$dispersion.scaled < Inf)
    top.features <- rownames(x = hvf.info)[which(x = means.use & dispersions.use)]
    VariableFeatures(object = assays[[assay]]) <- top.features
  } else if (verbose) {
    message("No variable feature expression found in h5ad file")
  }
  Key(object = assays[[assay]]) <- paste0(tolower(x = assay), '_')
  rm(x)
  CheckGC()
  # Get dimensional reduction information
  # If data isn't scaled, don't bother
  # if (scaled && file$exists(name = 'obsm')) {
  #   if (verbose) {
  #     message("Pulling dimensional reduction information")
  #     message("Pulling cell embeddings")
  #   }
  #   # Pull cell embeddings
  #   # embed.reduc <- file[['obsm']]$get_type()$get_cpd_labels()
  #   embed.reduc <- file[['obsm']] %>% names 
  #   embed.n <- sapply(
  #     X = embed.reduc, 
  #     FUN = function(tt) file[['obsm']][[tt]][,]
  #   )
  #   # embed.n <- sapply(
  #   #   X = file[['obsm']]$get_type()$describe()$cpd_types,
  #   #   FUN = '[[',
  #   #   'array_dims'
  #   # )
  #   names(x = embed.n) <- embed.reduc
  #   ncells <- file[['obsm']]$dims
  #   embeddings <- lapply(
  #     X = embed.reduc,
  #     FUN = function(r) {
  #       return(t(x = vapply(
  #         X = 1:ncells,
  #         FUN = function(i) {
  #           return(file[['obsm']][i][[r]])
  #         },
  #         FUN.VALUE = numeric(length = embed.n[[r]])
  #       )))
  #     }
  #   )
  #   names(x = embeddings) <- embed.reduc
  #   for (i in 1:length(x = embeddings)) {
  #     rownames(x = embeddings[[i]]) <- colnames(x = assays[[assay]])
  #   }
  #   # Pull feature loadings
  #   if (file$exists(name = 'varm')) {
  #     if (verbose) {
  #       message("Pulling feature loadings")
  #     }
  #     load.reduc <- file[['varm']]$get_type()$get_cpd_labels()
  #     load.n <- sapply(
  #       X = file[['varm']]$get_type()$describe()$cpd_types,
  #       FUN = '[[',
  #       'array_dims'
  #     )
  #     names(x = load.n) <- load.reduc
  #     nfeatures <- file[['varm']]$dims
  #     loadings <- lapply(
  #       X = load.reduc,
  #       FUN = function(r) {
  #         return(t(x = vapply(
  #           X = 1:nfeatures,
  #           FUN = function(i) {
  #             return(file[['varm']][i][[r]])
  #           },
  #           FUN.VALUE = numeric(length = load.n[[load.reduc]])
  #         )))
  #       }
  #     )
  #     match.ind <- lapply(
  #       X = gsub(pattern = 's$', replacement = '', x = tolower(x = load.reduc)),
  #       FUN = grep,
  #       x = embed.reduc
  #     )
  #     no.match <- which(x = sapply(X = match.ind, FUN = length) != 1)
  #     if (length(x = no.match) >= 1) {
  #       warning(
  #         "Unable to determine where the following feature loadings belong: ",
  #         paste(load.reduc[no.match], collapse = ', '),
  #         call. = FALSE,
  #         immediate. = TRUE
  #       )
  #       loadings <- loadings[-no.match]
  #       load.reduc <- load.reduc[-no.match]
  #       match.ind <- match.ind[-no.match]
  #     }
  #     names(x = loadings) <- embed.reduc[unlist(x = match.ind)]
  #     for (i in 1:length(x = loadings)) {
  #       rownames(x = loadings[[i]]) <- rownames(x = GetAssayData(
  #         object = assays[[assay]],
  #         slot = 'scale.data'
  #       ))
  #     }
  #   } else {
  #     if (verbose) {
  #       message("No feature loadings found")
  #     }
  #     loadings <- list()
  #   }
  #   # Create DimReduc objects
  #   dim.reducs <- vector(mode = 'list', length = length(x = embed.reduc))
  #   for (i in 1:length(x = embed.reduc)) {
  #     r <- embed.reduc[i]
  #     key <- tolower(x = gsub(pattern = 'X_', replacement = '', x = r))
  #     key <- switch(
  #       EXPR = key,
  #       'pca' = 'PC',
  #       'tsne' = 'tSNE',
  #       toupper(x = key)
  #     )
  #     key <- paste0(key, '_')
  #     stdev <- if (r == 'X_pca' && file$exists(name = 'uns') && file$exists(name = 'uns/pca/variance')) {
  #       sqrt(x = file[['uns/pca/variance']][])
  #     } else {
  #       numeric(length = 0L)
  #     }
  #     dim.reducs[[i]] <- CreateDimReducObject(
  #       embeddings = embeddings[[r]],
  #       loadings = loadings[[r]] %||% new(Class = 'matrix'),
  #       assay = assay,
  #       stdev = stdev,
  #       key = key
  #     )
  #   }
  #   # Properly name dimensional reductions
  #   names(x = dim.reducs) <- gsub(
  #     pattern = 'X_',
  #     replacement = '',
  #     x = embed.reduc
  #   )
  #   # Clean up
  #   rm(embeddings, loadings)
  #   CheckGC()
  # } else {
  #   if (verbose) {
  #     message("No dimensional reduction information found")
  #   }
  #   dim.reducs <- list()
  # }
  # # Create the Seurat object
  # if (verbose) {
  #   message("Assembling Seurat object")
  # }
  # Create a project name, will be used as identity classes
  project <- gsub(
    pattern = '\\.h5ad',
    replacement = '',
    x = basename(path = file$filename)
  )
  object <- new(
    Class = 'Seurat',
    assays = assays,
    meta.data = obs,
    version = packageVersion(pkg = 'Seurat'),
    project.name = project
  )
  # Set default assay and identity information
  DefaultAssay(object = object) <- assay
  Idents(object = object) <- project
  # # Add dimensional reduction infrom
  # if (scaled && length(x = dim.reducs) >= 1) {
  #   for (r in names(x = dim.reducs)) {
  #     object[[r]] <- dim.reducs[[r]]
  #   }
  # }
  # Get graph information
  # if (scaled && file$exists(name = 'uns') && file$exists(name = 'uns/neighbors')) {
  #   if (verbose) {
  #     message("Finding nearest neighbor graph")
  #   }
  #   graph <- as.sparse(x = file[['uns/neighbors/distances']])
  #   colnames(x = graph) <- rownames(x = graph) <- colnames(x = object)
  #   method <- ifelse(
  #     test = file[['uns/neighbors/params']]$exists(name = 'method'),
  #     yes = file[['uns/neighbors/params/method']][],
  #     no = 'adata'
  #   )
  #   object[[paste(assay, method, sep = '_')]] <- as.Graph(x = graph)
  # } else if (verbose) {
  #   message("No nearest-neighbor graph")
  # }
  return(object)
}
