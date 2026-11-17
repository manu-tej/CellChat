// Exact-parity sparse grouped Tukey trimean for CellChat permutation testing.
// Reproduces, bit for bit, mean(collapse::fquantile(x, c(.25,.5,.5,.75))) per (gene, group)
// where x is the group's slice of a gene row INCLUDING implicit zeros.
#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <cmath>
using namespace Rcpp;

// R's real_mean (summary.c) for the 4-vector (q25, q50, q50, q75), LDOUBLE = long double
static inline double r_mean4(double a, double b, double c, double d) {
  long double s = 0.0L;
  s += a; s += b; s += c; s += d;
  if (std::isfinite((double) s)) {
    s /= 4;
    if (std::isfinite((double) s)) {
      long double t = 0.0L;
      t += (a - s); t += (b - s); t += (c - s); t += (d - s);
      s += t / 4;
    }
  } else {
    s = 0.0L; s += a / 4.0; s += b / 4.0; s += c / 4.0; s += d / 4.0;  // overflow path (unused for [0,1] data)
    if (std::isfinite((double) s)) {
      long double t = 0.0L;
      t += (a - s) / 4; t += (b - s) / 4; t += (c - s) / 4; t += (d - s) / 4;
      s += t;
    }
  }
  return (double) s;
}

// Prepare: per gene (column of tX = t(data), a dgCMatrix N x G), nonzeros sorted ascending by value.
// [[Rcpp::export]]
List prep_sorted_nz(S4 tX) {
  IntegerVector p = tX.slot("p"), i = tX.slot("i");
  NumericVector x = tX.slot("x");
  int G = p.size() - 1;
  IntegerVector cell(x.size()); NumericVector val(x.size());
  std::vector<int> ord;
  for (int g = 0; g < G; ++g) {
    int st = p[g], en = p[g + 1], m = en - st;
    ord.resize(m);
    for (int k = 0; k < m; ++k) ord[k] = st + k;
    std::stable_sort(ord.begin(), ord.end(), [&](int a, int b) { return x[a] < x[b]; });
    for (int k = 0; k < m; ++k) { cell[st + k] = i[ord[k]]; val[st + k] = x[ord[k]]; }
  }
  return List::create(_["p"] = p, _["cell"] = cell, _["val"] = val, _["G"] = G);
}

// One grouping: grp is 0-based group per cell (length N), K groups. Returns G x K matrix of trimeans.
static void trimean_one(const int *p, const int *cell, const double *val, int G,
                        const int *grp, int K, const int *nk, double *out,
                        std::vector<int> &cnt, std::vector<int> &neg,
                        std::vector<int> &need, std::vector<double> &got, std::vector<int> &ptr,
                        std::vector<int> &nneed) {
  // per group, the 6 order-statistic ranks (ih, ih+1 for Q = .25, .5, .75) and fractions are fixed by n_k
  const double Qs[3] = {0.25, 0.5, 0.75};
  for (int g = 0; g < G; ++g) {
    int st = p[g], en = p[g + 1];
    std::fill(cnt.begin(), cnt.end(), 0);
    std::fill(neg.begin(), neg.end(), 0);
    for (int e = st; e < en; ++e) { int k = grp[cell[e]]; cnt[k]++; if (val[e] < 0) neg[k]++; }
    // needed nonzero ranks per group, in ascending order (map overall rank -> nonzero rank or zero)
    for (int k = 0; k < K; ++k) {
      int n = nk[k], z = n - cnt[k], ng = neg[k];
      nneed[k] = 0; ptr[k] = 0;
      for (int q = 0; q < 3; ++q) {
        double h = (n - 1) * Qs[q]; int ih = (int) h;
        for (int r = ih; r <= ih + 1 && r < n; ++r) {
          int nzr = (r < ng) ? r : ((r < ng + z) ? -1 : r - z);
          if (nzr >= 0) {
            int &m = nneed[k];
            if (m == 0 || need[k * 6 + m - 1] < nzr) { need[k * 6 + m] = nzr; m++; }
          }
        }
      }
    }
    std::fill(cnt.begin(), cnt.end(), 0);
    for (int e = st; e < en; ++e) {
      int k = grp[cell[e]]; int c = cnt[k]++;
      while (ptr[k] < nneed[k] && need[k * 6 + ptr[k]] == c) { got[k * 6 + ptr[k]] = val[e]; ptr[k]++; }
    }
    for (int k = 0; k < K; ++k) {
      int n = nk[k], z = n - cnt[k], ng = neg[k];
      auto ostat = [&](int r) -> double {
        if (r < ng) { for (int m = 0; m < nneed[k]; ++m) if (need[k * 6 + m] == r) return got[k * 6 + m]; }
        else if (r < ng + z) return 0.0;
        else { int nzr = r - z; for (int m = 0; m < nneed[k]; ++m) if (need[k * 6 + m] == nzr) return got[k * 6 + m]; }
        return NA_REAL;
      };
      double qv[3];
      for (int q = 0; q < 3; ++q) {
        if (n <= 1) { qv[q] = (n == 1) ? ostat(0) : NA_REAL; continue; }
        double h = (n - 1) * Qs[q]; int ih = (int) h; double fr = h - ih;   // collapse RETQSWITCH case 7
        double a = ostat(ih);
        if (ih == n - 1 || fr <= 0.0) qv[q] = a;                               // dquickselect_elem
        else { double b = ostat(ih + 1); qv[q] = a + fr * (b - a); }
      }
      out[g + (size_t) k * G] = r_mean4(qv[0], qv[1], qv[1], qv[2]);
    }
  }
}

// perm: N x B integer matrix of 1-based cell indices (CellChat's `permutation`); group: 1-based ints of length N.
// Returns a list of B matrices (G x K), the same as data.use.avg.boot in computeCommunProb.
// [[Rcpp::export]]
List trimean_perm(List prep, IntegerVector group, int K, IntegerMatrix perm, bool identity = false) {
  IntegerVector p = prep["p"], cell = prep["cell"]; NumericVector val = prep["val"]; int G = prep["G"];
  int N = group.size();
  std::vector<int> nk(K, 0); for (int c = 0; c < N; ++c) nk[group[c] - 1]++;
  std::vector<int> cnt(K), neg(K), need(6 * K), ptr(K), nneed(K); std::vector<double> got(6 * K);
  std::vector<int> grp(N);
  int B = identity ? 1 : perm.ncol();
  List res(B);
  for (int b = 0; b < B; ++b) {
    if (identity) for (int c = 0; c < N; ++c) grp[c] = group[c] - 1;
    else for (int c = 0; c < N; ++c) grp[c] = group[perm(c, b) - 1] - 1;   // groupboot <- group[permutation[, nE]]
    NumericMatrix out(G, K);
    trimean_one(p.begin(), cell.begin(), val.begin(), G, grp.data(), K, nk.data(), out.begin(), cnt, neg, need, got, ptr, nneed);
    res[b] = out;
  }
  return res;
}
