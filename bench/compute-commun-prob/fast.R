computeCommunProb_fast <- function(object, type = c("triMean", "truncatedMean","thresholdedMean", "median"), trim = 0.1, LR.use = NULL, raw.use = TRUE, population.size = FALSE,
                              distance.use = TRUE, interaction.range = 250, scale.distance = 0.01, k.min = 10, contact.dependent = TRUE, contact.range = NULL, contact.knn.k = NULL, contact.dependent.forced = FALSE, do.symmetric = TRUE,
                              nboot = 100, seed.use = 1L, Kh = 0.5, n = 1) {
  type <- match.arg(type)
  cat(type, "is used for calculating the average gene expression per cell group.", "\n")
  FunMean <- switch(type,
                    triMean = triMean,
                    truncatedMean = function(x) mean(x, trim = trim, na.rm = TRUE),
                    thresholdedMean = function(x) thresholdedMean(x, trim = trim, na.rm = TRUE),
                    median = function(x) median(x, na.rm = TRUE))

  fast.trimean <- (type == "triMean")
  if (raw.use) {
    data <- object@data.signaling
  } else {
    data <- object@data.smooth
  }
  if (fast.trimean) {
    data <- methods::as(methods::as(methods::as(data, "dMatrix"), "generalMatrix"), "CsparseMatrix")
    if (anyNA(data@x)) { fast.trimean <- FALSE; data <- as.matrix(data) }
  } else {
    data <- as.matrix(data)
  }
  if (is.null(LR.use)) {
    pairLR.use <- object@LR$LRsig
  } else {
    if (length(unique(LR.use$annotation)) > 1) {
      LR.use$annotation <- factor(LR.use$annotation, levels = c("Secreted Signaling","ECM-Receptor", "Non-protein Signaling", "Cell-Cell Contact"))
      LR.use <- LR.use[order(LR.use$annotation), , drop = FALSE]
      LR.use$annotation <- as.character(LR.use$annotation)
    }
    pairLR.use <- LR.use
  }
  complex_input <- object@DB$complex
  cofactor_input <- object@DB$cofactor
  my.sapply <- ifelse(
    test = future::nbrOfWorkers() == 1,
    yes = sapply,
    no = future.apply::future_sapply
  )

  ptm = Sys.time()

  pairLRsig <- pairLR.use
  group <- object@idents
  geneL <- as.character(pairLRsig$ligand)
  geneR <- as.character(pairLRsig$receptor)
  nLR <- nrow(pairLRsig)
  numCluster <- nlevels(group)
  if (numCluster != length(unique(group))) {
    stop("Please check `unique(object@idents)` and ensure that the factor levels are correct!
         You may need to drop unused levels using 'droplevels' function. e.g.,
         `meta$labels = droplevels(meta$labels, exclude = setdiff(levels(meta$labels),unique(meta$labels)))`")
  }

  if (fast.trimean) {
    data.use <- data; data.use@x <- data@x / max(data)   # same values as dense data/max(data)
    data.use <- Matrix::drop0(data.use)
  } else {
    data.use <- data/max(data)
  }
  nC <- ncol(data.use)

  # compute the average expression per group
  if (fast.trimean) {
    prep.tm <- prep_sorted_nz(methods::as(Matrix::t(data.use), "CsparseMatrix"))
    data.use.avg <- trimean_perm(prep.tm, as.integer(group), numCluster, matrix(0L, 1, 1), identity = TRUE)[[1]]
    rownames(data.use.avg) <- rownames(data.use)
  } else {
    data.use.avg <- aggregate(t(data.use), list(group), FUN = FunMean)
    data.use.avg <- t(data.use.avg[,-1])
  }
  colnames(data.use.avg) <- levels(group)
  # compute the expression of ligand or receptor
  dataLavg <- computeExpr_LR(geneL, data.use.avg, complex_input)
  dataRavg <- computeExpr_LR(geneR, data.use.avg, complex_input)
  # take account into the effect of co-activation and co-inhibition receptors
  dataRavg.co.A.receptor <- computeExpr_coreceptor(cofactor_input, data.use.avg, pairLRsig, type = "A")
  dataRavg.co.I.receptor <- computeExpr_coreceptor(cofactor_input, data.use.avg, pairLRsig, type = "I")
  dataRavg <- dataRavg * dataRavg.co.A.receptor/dataRavg.co.I.receptor

  dataLavg2 <- t(replicate(nrow(dataLavg), as.numeric(table(group))/nC))
  dataRavg2 <- dataLavg2

  # compute the expression of agonist and antagonist
  index.agonist <- which(!is.na(pairLRsig$agonist) & pairLRsig$agonist != "")
  index.antagonist <- which(!is.na(pairLRsig$antagonist) & pairLRsig$antagonist != "")
  # quantify the communication probability

  # compute the spatial constraint
  if (object@options$datatype != "RNA") {
    data.spatial <- object@images$coordinates
    if ("spatial.factors" %in% names(object@images)) {
      ratio <- object@images$spatial.factors$ratio
      tol <- object@images$spatial.factors$tol
    } else {
      stop("`object@images$spatial.factors` is missing. Please update the object via `updateCellChat`! \n")
    }

    meta.t = data.frame(group = group, samples = object@meta$samples, row.names = rownames(object@meta))
    res <- computeRegionDistance(coordinates = data.spatial, meta = meta.t, interaction.range = interaction.range, ratio = ratio, tol = tol, k.min = k.min, contact.dependent = contact.dependent, contact.range = contact.range, contact.knn.k = contact.knn.k)
    d.spatial <- res$d.spatial # NaN if no nearby cell pairs
    adj.contact <- res$adj.contact # zeros if no nearby cell pairs
    if (distance.use) {
      print(paste0('>>> Run CellChat on spatial transcriptomics data using distances as constraints of the computed communication probability <<< [', Sys.time(),']'))
      d.spatial <- d.spatial * scale.distance
      diag(d.spatial) <- NaN
      d.min <- min(d.spatial, na.rm = TRUE)
      if (d.min < 1) {
        cat("The suggested minimum value of scaled distances is in [1,2], and the calculated value here is ", d.min,"\n")
        stop("Please increase the value of `scale.distance` and use a value that is slighly smaller than ", format(1/d.min, digits = 2) ,"\n")
      }
      P.spatial <- 1/d.spatial
      P.spatial[is.na(d.spatial)] <- 0
      diag(P.spatial) <- max(P.spatial) # if this value is 1, the self-connections will have more larger weight.
      d.spatial <- d.spatial/scale.distance # This is only for saving the data
    } else {
      print(paste0('>>> Run CellChat on spatial transcriptomics data without distance values as constraints of the computed communication probability <<< [', Sys.time(),']'))
      P.spatial <- matrix(1, nrow = numCluster, ncol = numCluster)
      P.spatial[is.na(d.spatial)] <- 0 # diagonal is 1
    }

  } else {
    print(paste0('>>> Run CellChat on sc/snRNA-seq data <<< [', Sys.time(),']'))
    d.spatial <- matrix(NaN, nrow = numCluster, ncol = numCluster)
    P.spatial <- matrix(1, nrow = numCluster, ncol = numCluster)
    adj.contact <- matrix(1, nrow = numCluster, ncol = numCluster)
    contact.dependent = FALSE; contact.dependent.forced = FALSE; contact.range = NULL; contact.knn.k = NULL;
    distance.use = NULL; interaction.range = NULL; ratio = NULL; tol = NULL; k.min = NULL;
  }

  if (object@options$datatype == "RNA") {
    nLR1 <- nLR
  } else {
    if (contact.dependent.forced == TRUE) {
      cat("Force to run CellChat in a `contact-dependent` manner for all L-R pairs including secreted signaling.\n")
      P.spatial <- P.spatial * adj.contact
      nLR1 <- nLR
    } else { # contact.dependent.forced == F
      if (contact.dependent == TRUE && length(unique(pairLRsig$annotation)) > 0) {
        if (all(unique(pairLRsig$annotation) %in% c("Cell-Cell Contact"))) {
          cat("All the input L-R pairs are `Cell-Cell Contact` signaling. Run CellChat in a contact-dependent manner. \n")
          P.spatial <- P.spatial * adj.contact
          nLR1 <- nLR
        } else if (all(unique(pairLRsig$annotation) %in% c("Secreted Signaling", "ECM-Receptor", "Non-protein Signaling"))) {
          cat("Molecules of the input L-R pairs are diffusible. Run CellChat in a diffusion manner based on the `interaction.range`.\n")
          nLR1 <- nLR
        } else {
          cat("The input L-R pairs have both secreted signaling and contact-dependent signaling. Run CellChat in a contact-dependent manner for `Cell-Cell Contact` signaling, and in a diffusion manner based on the `interaction.range` for other L-R pairs. \n")
          nLR1 <- max(which(pairLRsig$annotation %in% c("Secreted Signaling", "ECM-Receptor", "Non-protein Signaling")))
        }
      } else { # contact.dependent == F or there is no `annotation` column in the database
        cat("Run CellChat in a diffusion manner based on the `interaction.range` for all L-R pairs. Setting `contact.dependent = TRUE` if preferring a contact-dependent manner for `Cell-Cell Contact` signaling. \n")
        nLR1 <- nLR
      }
    }
  }

  Prob <- array(0, dim = c(numCluster,numCluster,nLR))
  Pval <- array(0, dim = c(numCluster,numCluster,nLR))

  set.seed(seed.use)
  permutation <- replicate(nboot, sample.int(nC, size = nC))
  if (fast.trimean) {
    data.use.avg.boot <- trimean_perm(prep.tm, as.integer(group), numCluster, permutation)
    data.use.avg.boot <- lapply(data.use.avg.boot, function(m) { rownames(m) <- rownames(data.use); m })
  } else {
  data.use.avg.boot <- my.sapply(
    X = 1:nboot,
    FUN = function(nE) {
      groupboot <- group[permutation[, nE]]
      data.use.avgB <- aggregate(t(data.use), list(groupboot), FUN = FunMean)
      data.use.avgB <- t(data.use.avgB[,-1])
      return(data.use.avgB)
    },
    simplify = FALSE
  )
  }
  # ---- permutations batched: every L-R expression is a per-column function of the group-average matrix,
  # so evaluate each unique ligand/receptor/cofactor ONCE on the G x (K*nboot) matrix of all permuted averages ----
  AvgAll <- do.call(cbind, data.use.avg.boot)          # column block nE = permutation nE
  uL <- unique(geneL); uR <- unique(geneR)
  EL <- computeExpr_LR(uL, AvgAll, complex_input); ER <- computeExpr_LR(uR, AvgAll, complex_input)
  jL <- match(geneL, uL); jR <- match(geneR, uR)
  uA <- unique(pairLRsig$co_A_receptor); uI <- unique(pairLRsig$co_I_receptor)
  CA <- computeExpr_coreceptor(cofactor_input, AvgAll, data.frame(co_A_receptor = uA, stringsAsFactors = FALSE), type = "A")
  CI <- computeExpr_coreceptor(cofactor_input, AvgAll, data.frame(co_I_receptor = uI, stringsAsFactors = FALSE), type = "I")
  jA <- match(pairLRsig$co_A_receptor, uA); jI <- match(pairLRsig$co_I_receptor, uI)
  ag.names <- unique(pairLRsig$agonist[index.agonist]); an.names <- unique(pairLRsig$antagonist[index.antagonist])
  ag.first <- index.agonist[match(ag.names, pairLRsig$agonist[index.agonist])]
  an.first <- index.antagonist[match(an.names, pairLRsig$antagonist[index.antagonist])]
  AG <- lapply(ag.first, function(a) as.vector(computeExpr_agonist(data.use = AvgAll, pairLRsig, cofactor_input, index.agonist = a, Kh = Kh, n = n)))
  AN <- lapply(an.first, function(a) as.vector(computeExpr_antagonist(data.use = AvgAll, pairLRsig, cofactor_input, index.antagonist = a, Kh = Kh, n = n)))
  rm(AvgAll)
  ik <- rep(1:numCluster, times = numCluster); il <- rep(1:numCluster, each = numCluster)
  P4.vec <- if (population.size) as.vector(Matrix::crossprod(matrix(as.numeric(table(group))/nC, nrow = 1))) else NULL
  pb <- txtProgressBar(min = 0, max = nLR, style = 3, file = stderr())

  for (i in 1:nLR) {
    # ligand/receptor
    dataLR <- Matrix::crossprod(matrix(dataLavg[i,], nrow = 1), matrix(dataRavg[i,], nrow = 1))
    P1 <- dataLR^n/(Kh^n + dataLR^n)
    P1_Pspatial <- P1*P.spatial
    if (sum(P1_Pspatial) == 0) {
      Pnull = P1_Pspatial
      Prob[ , , i] <- Pnull
      p = 1
      Pval[, , i] <- matrix(p, nrow = numCluster, ncol = numCluster, byrow = FALSE)
    } else {
      if (i > nLR1) {
        P.spatial <- P.spatial * adj.contact
      }
      # agonist and antagonist
      if (is.element(i, index.agonist)) {
        data.agonist <- computeExpr_agonist(data.use = data.use.avg, pairLRsig, cofactor_input, index.agonist = i, Kh = Kh,  n = n)
        P2 <- Matrix::crossprod(matrix(data.agonist, nrow = 1))
      } else {
        P2 <- matrix(1, nrow = numCluster, ncol = numCluster)
      }
      if (is.element(i, index.antagonist)) {
        data.antagonist <- computeExpr_antagonist(data.use = data.use.avg, pairLRsig, cofactor_input,  index.antagonist = i, Kh = Kh,  n = n)
        P3 <- Matrix::crossprod(matrix(data.antagonist, nrow = 1))
      } else {
        P3 <- matrix(1, nrow = numCluster, ncol = numCluster)
      }
      # number of cells
      if (population.size) {
        P4 <- Matrix::crossprod(matrix(dataLavg2[i,], nrow = 1), matrix(dataRavg2[i,], nrow = 1))
      } else {
        P4 <- matrix(1, nrow = numCluster, ncol = numCluster)
      }

      # Pnull = P1*P2*P3*P4
      Pnull = P1*P2*P3*P4*P.spatial
      Prob[ , , i] <- Pnull

      Pnull <- as.vector(Pnull)

      Lb <- matrix(EL[jL[i], ], nrow = numCluster)
      Rb <- matrix(ER[jR[i], ] * CA[jA[i], ]/CI[jI[i], ], nrow = numCluster)
      dataLRB <- Lb[ik, , drop = FALSE] * Rb[il, , drop = FALSE]
      Pboot <- dataLRB^n/(Kh^n + dataLRB^n)
      if (is.element(i, index.agonist)) {
        A <- matrix(AG[[match(pairLRsig$agonist[i], ag.names)]], nrow = numCluster)
        Pboot <- Pboot * (A[ik, , drop = FALSE] * A[il, , drop = FALSE])
      }
      if (is.element(i, index.antagonist)) {
        A <- matrix(AN[[match(pairLRsig$antagonist[i], an.names)]], nrow = numCluster)
        Pboot <- Pboot * (A[ik, , drop = FALSE] * A[il, , drop = FALSE])
      }
      if (population.size) Pboot <- Pboot * P4.vec
      Pboot <- Pboot * as.vector(P.spatial)
      nReject <- rowSums(Pboot - Pnull > 0)
      p = nReject/nboot
      Pval[, , i] <- matrix(p, nrow = numCluster, ncol = numCluster, byrow = FALSE)
    }
    setTxtProgressBar(pb = pb, value = i)
  }
  close(con = pb)
  Pval[Prob == 0] <- 1
  dimnames(Prob) <- list(levels(group), levels(group), rownames(pairLRsig))
  dimnames(Pval) <- dimnames(Prob)
  net <- list("prob" = Prob, "pval" = Pval)
  execution.time = Sys.time() - ptm
  object@options$run.time <- as.numeric(execution.time, units = "secs")

  object@options$parameter <- list(type.mean = type, trim = trim, raw.use = raw.use, population.size = population.size,  nboot = nboot, seed.use = seed.use, Kh = Kh, n = n,
                                   distance.use = distance.use, interaction.range = interaction.range, ratio = ratio, tol = tol, k.min = k.min,
                                   contact.dependent = contact.dependent, contact.range = contact.range, contact.knn.k = contact.knn.k, contact.dependent.forced = contact.dependent.forced
                                   )
  if (object@options$datatype != "RNA") {
    object@images$distance <- d.spatial
  }
  object@net <- net
  print(paste0('>>> CellChat inference is done. Parameter values are stored in `object@options$parameter` <<< [', Sys.time(),']'))
  return(object)
}
