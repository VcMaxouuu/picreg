// Harrell's concordance index in O(n log n).
//
// A pair (i, j) is comparable when i is an event (events[i] == 1) and either
// t_i < t_j, or t_i == t_j and j is censored (events[j] == 0). It is
// concordant when pred_i > pred_j and counts 1/2 when pred_i == pred_j.
//
// Subjects are swept by decreasing time. A Fenwick tree over the ranks of
// the predictions holds every subject with a strictly later time; the
// censored subjects of the current tie group are added before the events of
// that group are queried, and the remaining ones after.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include <algorithm>
#include <numeric>
#include <vector>

// [[Rcpp::export]]
double concordance_index_cpp(const arma::vec& times,
                             const arma::vec& events,
                             const arma::vec& predictions) {
  const std::size_t n = times.n_elem;

  // Dense 1-based ranks of the predictions; ties share a rank.
  std::vector<std::size_t> by_pred(n);
  std::iota(by_pred.begin(), by_pred.end(), 0);
  std::sort(by_pred.begin(), by_pred.end(),
            [&predictions](std::size_t a, std::size_t b) {
              return predictions[a] < predictions[b];
            });
  std::vector<std::size_t> rank(n);
  std::size_t n_ranks = 0;
  for (std::size_t k = 0; k < n; ++k) {
    if (k == 0 || predictions[by_pred[k]] != predictions[by_pred[k - 1]])
      ++n_ranks;
    rank[by_pred[k]] = n_ranks;
  }

  // Fenwick tree of subject counts per prediction rank.
  std::vector<double> tree(n_ranks + 1, 0.0);
  auto add = [&tree, n_ranks](std::size_t i) {
    for (; i <= n_ranks; i += i & (~i + 1)) tree[i] += 1.0;
  };
  auto count_up_to = [&tree](std::size_t i) {
    double s = 0.0;
    for (; i > 0; i -= i & (~i + 1)) s += tree[i];
    return s;
  };

  std::vector<std::size_t> by_time(n);
  std::iota(by_time.begin(), by_time.end(), 0);
  std::sort(by_time.begin(), by_time.end(),
            [&times](std::size_t a, std::size_t b) {
              return times[a] > times[b];
            });

  double concordant = 0.0, tied = 0.0, comparable = 0.0, in_tree = 0.0;
  std::size_t start = 0;
  while (start < n) {
    std::size_t end = start;
    while (end < n && times[by_time[end]] == times[by_time[start]]) ++end;

    for (std::size_t k = start; k < end; ++k) {
      if (events[by_time[k]] == 0.0) {
        add(rank[by_time[k]]);
        in_tree += 1.0;
      }
    }
    for (std::size_t k = start; k < end; ++k) {
      const std::size_t i = by_time[k];
      if (events[i] == 1.0) {
        const double below = count_up_to(rank[i] - 1);
        concordant += below;
        tied       += count_up_to(rank[i]) - below;
        comparable += in_tree;
      }
    }
    for (std::size_t k = start; k < end; ++k) {
      if (events[by_time[k]] != 0.0) {
        add(rank[by_time[k]]);
        in_tree += 1.0;
      }
    }
    start = end;
  }

  if (comparable == 0.0) return 0.0;
  return (concordant + 0.5 * tied) / comparable;
}
