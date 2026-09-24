// TestJetSorter
// pT order maps on hand-made arrays, ties, reordering, and random arrays

#include "JetSorter.h"

#include <cstdio>
#include <random>
#include <string>
#include <vector>

static int nFail = 0;
static int nCheck = 0;

static void Check(bool ok, const std::string &what) {
  nCheck++;
  if (!ok) {
    nFail++;
    std::printf("  FAIL: %s\n", what.c_str());
  }
}

static void TestOrder() {
  std::printf("order map\n");
  const float pt[4] = {50, 80, 40, 30};
  const float eta[4] = {0.1f, -1.2f, 2.0f, 0.5f};
  std::vector<int> order = JetSorter::Order(4, pt);
  Check(order == std::vector<int>({1, 0, 2, 3}), "{50,80,40,30} -> {1,0,2,3}");
  Check(eta[order[0]] == -1.2f, "leading eta read through the map");

  std::vector<float> etaSorted = JetSorter::Reorder(order, eta);
  Check(etaSorted == std::vector<float>({-1.2f, 0.1f, 2.0f, 0.5f}),
        "Reorder copies eta into pT order");

  const double ties[4] = {10, 20, 20, 5};
  Check(JetSorter::Order(4, ties) == std::vector<int>({1, 2, 0, 3}),
        "equal pT keep their original order (double input)");
  Check(JetSorter::Order(0, pt).empty(), "no jets, empty map");
}

static void TestRandom() {
  std::printf("random arrays\n");
  std::mt19937 rng(7);
  std::uniform_real_distribution<float> u(0, 500);
  bool ok = true;
  for (int trial = 0; trial < 1000; trial++) {
    const int n = trial % 40;
    std::vector<float> pt(n);
    for (auto &p : pt) {
      p = u(rng);
    }
    std::vector<int> order = JetSorter::Order(n, pt.data());
    std::vector<bool> seen(n, false);
    for (int i = 0; i < n; i++) {
      ok = ok && order[i] >= 0 && order[i] < n && !seen[order[i]];
      seen[order[i]] = true;
      ok = ok && (i == 0 || pt[order[i - 1]] >= pt[order[i]]);
    }
  }
  Check(ok, "1000 random arrays: permutation, pT descending");
}

int main() {
  TestOrder();
  TestRandom();
  std::printf("%d/%d checks passed\n", nCheck - nFail, nCheck);
  return nFail == 0 ? 0 : 1;
}
