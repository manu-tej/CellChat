make_input <- function(N, K, nLR = NULL, seed = 42, density_base = c(0.5, 10)) {
  set.seed(seed)
  db <- CellChatDB.human
  inter <- db$interaction; cpx <- db$complex; cof <- db$cofactor
  cpx_genes <- unlist(cpx[, grepl("subunit", colnames(cpx))]); cpx_genes <- cpx_genes[cpx_genes != ""]
  lr_single <- setdiff(unique(c(inter$ligand, inter$receptor)), rownames(cpx))
  cof_genes <- unlist(cof[, grepl("cofactor", colnames(cof))]); cof_genes <- cof_genes[cof_genes != ""]
  genes <- sort(unique(c(lr_single, cpx_genes, cof_genes)))
  G <- length(genes)
  # group sizes: imbalanced, like a tumour atlas
  w <- rgamma(K, shape = 1.2); w <- w / sum(w)
  sizes <- pmax(20, round(w * N)); sizes[1] <- sizes[1] + (N - sum(sizes))
  grp <- factor(rep(paste0("C", seq_len(K)), sizes), levels = paste0("C", seq_len(K)))
  grp <- grp[sample.int(N)]
  # detection rates per gene x group: low baseline, a few "marker" groups with high detection
  p <- matrix(rbeta(G * K, density_base[1], density_base[2]), G, K)
  nmark <- rpois(G, 1.5)
  for (g in seq_len(G)) if (nmark[g] > 0) p[g, sample.int(K, min(K, nmark[g]))] <- runif(min(K, nmark[g]), 0.2, 0.9)
  gi <- as.integer(grp)
  # simulate sparse matrix column by column chunks
  ii <- list(); jj <- list(); xx <- list()
  for (k in seq_len(K)) {
    cells <- which(gi == k); nk <- length(cells)
    for (g in seq_len(G)) {
      m <- rbinom(1, nk, p[g, k]); if (m == 0) next
      cc <- cells[sample.int(nk, m)]
      ii[[length(ii) + 1]] <- rep.int(g, m); jj[[length(jj) + 1]] <- cc
      xx[[length(xx) + 1]] <- log1p(rgamma(m, shape = 1.5, rate = 0.6))
    }
  }
  X <- sparseMatrix(i = unlist(ii), j = unlist(jj), x = unlist(xx), dims = c(G, N), dimnames = list(genes, paste0("cell", seq_len(N))))
  # "over-expressed" genes: detection > 0.2 in some group (stand-in for identifyOverExpressedGenes)
  feat <- genes[apply(p, 1, max) > 0.2]
  cpx_ok <- rownames(cpx)[sapply(seq_len(nrow(cpx)), function(r) { s <- unlist(cpx[r, ]); s <- s[s != ""]; any(s %in% feat) && all(s %in% genes) })]
  keep <- (inter$ligand %in% c(feat, cpx_ok)) & (inter$receptor %in% c(feat, cpx_ok))
  LRsig <- inter[keep, ]
  if (!is.null(nLR) && nLR < nrow(LRsig)) LRsig <- LRsig[sort(sample.int(nrow(LRsig), nLR)), ]
  obj <- new("CellChatLite", data.signaling = X, data.smooth = X, LR = list(LRsig = LRsig), DB = db,
             idents = grp, options = list(datatype = "RNA"), images = list(),
             meta = data.frame(row.names = colnames(X)), net = list())
  obj
}
