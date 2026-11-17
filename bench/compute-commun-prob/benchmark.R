args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L || !args[[1L]] %in% c("prepare", "main", "prototype", "compare") ||
    !args[[2L]] %in% c("original", "larger")) {
  stop("Usage: Rscript bench/compute-commun-prob/benchmark.R MODE CASE OUTPUT_DIR",
       "\nMODE: prepare, main, prototype, or compare. CASE: original or larger.", call. = FALSE)
}
mode <- args[[1L]]
case <- args[[2L]]
dir.create(args[[3L]], recursive = TRUE, showWarnings = FALSE)
results <- normalizePath(args[[3L]], mustWork = TRUE)
input <- file.path(results, paste0(case, "-input.rds"))
code <- file.path("bench", "compute-commun-prob")

if (mode == "compare") {
  reference <- readRDS(file.path(results, paste0(case, "-main.rds")))
  prototype <- readRDS(file.path(results, paste0(case, "-prototype.rds")))
  for (metric in c("prob", "pval")) {
    a <- reference[[metric]]
    b <- prototype[[metric]]
    stopifnot(identical(dim(a), dim(b)), identical(dimnames(a), dimnames(b)),
              all(is.finite(a)), all(is.finite(b)))
    delta <- abs(a - b)
    cat(sprintf("%s %s values=%d identical=%s max_abs_delta=%.17g nonzero_deltas=%d\n",
                case, metric, length(a), identical(a, b), max(delta), sum(delta != 0)))
  }
  stopifnot(identical(reference$prob, prototype$prob), identical(reference$pval, prototype$pval))
  quit(save = "no")
}

suppressPackageStartupMessages(library(Matrix))
setClass("CellChatLite", representation(data.signaling = "ANY", data.smooth = "ANY", LR = "list",
         DB = "list", idents = "factor", options = "list", images = "list", meta = "data.frame", net = "list"))

if (mode == "prepare") {
  load("data/CellChatDB.human.rda")
  source(file.path(code, "make-input.R"))
  object <- switch(case, original = make_input(20000L, 20L), larger = make_input(60000L, 65L, 3055L))
  if (case == "original") stopifnot(nrow(object@LR$LRsig) == 2719L)
  if (case == "larger") stopifnot(nrow(object@LR$LRsig) == 3055L)
  cat(sprintf("%s: %d cells, %d genes, %d groups, %d pairs, %d nonzero expression values\n",
              case, ncol(object@data.signaling), nrow(object@data.signaling), nlevels(object@idents),
              nrow(object@LR$LRsig), Matrix::nnzero(object@data.signaling)))
  saveRDS(object, input, compress = FALSE)
  quit(save = "no")
}

suppressPackageStartupMessages({library(dplyr); library(magrittr); library(future)})
source("R/modeling.R")
future::plan(future::sequential)
run <- computeCommunProb
if (mode == "prototype") {
  cache <- file.path(results, "cppcache")
  dir.create(cache, showWarnings = FALSE)
  Rcpp::sourceCpp(file.path(code, "grouped_trimean.cpp"), cacheDir = cache)
  source(file.path(code, "fast.R"))
  run <- computeCommunProb_fast
}
object <- readRDS(input)
cat(sprintf("Running %s %s at %s\n", case, mode, format(Sys.time())))
flush.console()
gc()
elapsed <- system.time({ output <- run(object, type = "triMean", nboot = 100L, seed.use = 1L) })[["elapsed"]]
cat(sprintf("RESULT %s %s elapsed_s=%.3f nonzero_pairs=%d at %s\n", case, mode, elapsed,
            sum(colSums(matrix(output@net$prob, ncol = nrow(object@LR$LRsig))) > 0), format(Sys.time())))
saveRDS(output@net, file.path(results, paste0(case, "-", mode, ".rds")), compress = FALSE)
