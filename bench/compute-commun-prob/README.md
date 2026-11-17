# `computeCommunProb()` benchmark prototype

This standalone harness sources `R/modeling.R` from the checkout and uses the bundled CellChat database to generate two synthetic inputs. The prototype uses sparse grouped trimeans and reuses ligand, receptor, and cofactor calculations across pairs in the bootstrap. The baseline is CellChat `main` at `75253cd0c9e68410e6e721a6d3a0419a1d7e358f`.

Install `Matrix`, `dplyr`, `magrittr`, `future`, `Rcpp`, and `collapse` in R and have a C++ toolchain available for `Rcpp::sourceCpp()`. The measured runs used R 4.5.1, Matrix 1.7-3, dplyr 1.2.1, future 1.82.0, Rcpp 1.1.2, and collapse 2.1.9.

From the repository root, run both cases:

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
results_dir=$(mktemp -d)
for benchmark_case in original larger; do
  Rscript --vanilla bench/compute-commun-prob/benchmark.R prepare "$benchmark_case" "$results_dir" || exit
  Rscript --vanilla bench/compute-commun-prob/benchmark.R main "$benchmark_case" "$results_dir" || exit
  Rscript --vanilla bench/compute-commun-prob/benchmark.R prototype "$benchmark_case" "$results_dir" || exit
  Rscript --vanilla bench/compute-commun-prob/benchmark.R compare "$benchmark_case" "$results_dir" || exit
done
```

The two cases use 20,000 cells across 20 groups with 2,719 pairs, and 60,000 cells across 65 groups with 3,055 pairs. Each timed run uses one worker, `triMean`, 100 bootstraps, and `seed.use = 1`. Timing covers `computeCommunProb()` or its prototype. The comparison reports the maximum absolute difference for each output and checks exact equality of the probability and p-value arrays.
