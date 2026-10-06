// Two-dimensional binning for adaptive broadening (c++/adaptive-bins.hpp).

#include <cmath>
#include <complex>
#include <cstdio>
#include <fstream>
#include <limits>
#include <string>
#include <vector>
#include <gtest/gtest.h>

#include <params.hpp>
#include <bins.hpp>
#include <adaptive-bins.hpp>

using namespace NRG;

namespace {

void setup(Params &P) {
  P.emin = 1e-6;
  P.emax = 10.0;
  P.bins = 10;
  P.adaptive_broadening_sigma_step = 0.01;
  P.adaptive_broadening_sigma_max = 1.0;
}

const std::vector<double> energies = {3e-7, 1e-6, 2.5e-6, 1.234e-3, 0.5, 0.999, 9.9, 50.0};

} // namespace

TEST(AdaptiveBins, grid_matches_Bins) { // NOLINT
  Params P;
  setup(P);
  Bins<double> b(P);
  AdaptiveBins<double> a(P);
  ASSERT_EQ(a.omega_points(), b.bins.size());
  for (std::size_t i = 0; i < a.omega_points(); i++)
    EXPECT_EQ(a.omega(i), b.bins[i].first);
  EXPECT_EQ(a.sigma_points(), 100);
  EXPECT_DOUBLE_EQ(a.sigma(0), 0.01);
  EXPECT_DOUBLE_EQ(a.sigma(99), 1.0);
}

TEST(AdaptiveBins, omega_split_matches_Bins) { // NOLINT
  Params P;
  setup(P);
  Bins<double> b(P);
  AdaptiveBins<double> a(P);
  for (std::size_t n = 0; n < energies.size(); n++) {
    const std::complex<double> w(1.0 + 0.1 * double(n), 0.0);
    b.add(energies[n], w);
    a.add(energies[n], 0.05 * double(n + 1), w); // a different sigma per weight
  }
  // Summed over sigma, the 2D bins equal the 1D bins.
  std::vector<std::complex<double>> summed(a.omega_points());
  for (const auto &c : a.cells()) {
    const auto i = std::size_t(std::llround(std::log10(c.energy / a.omega(0)) * P.bins));
    ASSERT_LT(i, summed.size());
    EXPECT_DOUBLE_EQ(c.energy, a.omega(i));
    summed[i] += c.weight;
  }
  for (std::size_t i = 0; i < summed.size(); i++)
    EXPECT_NEAR(std::abs(summed[i] - b.bins[i].second), 0.0, 1e-14) << "omega point " << i;
  EXPECT_NEAR(std::abs(a.total_weight() - b.total_weight()), 0.0, 1e-13);
}

TEST(AdaptiveBins, sigma_nearest_and_clamped) { // NOLINT
  Params P;
  setup(P);
  AdaptiveBins<double> a(P);
  EXPECT_EQ(a.sigma_index(0.01), 0);
  EXPECT_EQ(a.sigma_index(0.0149), 0);
  EXPECT_EQ(a.sigma_index(0.0151), 1);
  EXPECT_EQ(a.sigma_index(0.693), 68);
  EXPECT_EQ(a.sigma_index(0.0), 0);   // below the grid
  EXPECT_EQ(a.sigma_index(-1.0), 0);
  EXPECT_EQ(a.sigma_index(1.0), 99);  // the last point
  EXPECT_EQ(a.sigma_index(7.5), 99);  // beyond the grid
  EXPECT_EQ(a.sigma_index(std::numeric_limits<double>::infinity()), 99);
  EXPECT_EQ(a.sigma_index(std::numeric_limits<double>::quiet_NaN()), 99);
  a.add(0.5, 0.6931, 1.0);
  const auto cells = a.cells();
  ASSERT_FALSE(cells.empty());
  for (const auto &c : cells) EXPECT_DOUBLE_EQ(c.sigma, 0.69);
}

TEST(AdaptiveBins, discard_immediately) { // NOLINT
  Params P;
  setup(P);
  P.discard_immediately = 1e-3;
  AdaptiveBins<double> a(P);
  a.add(1.0, 0.5, 1e-4); // below 1e-3 * energy
  EXPECT_EQ(a.nr_cells(), 0);
  a.add(1.0, 0.5, 1e-2);
  EXPECT_GT(a.nr_cells(), 0);
}

TEST(AdaptiveBins, merge) { // NOLINT
  Params P;
  setup(P);
  AdaptiveBins<double> a(P), b(P);
  a.add(0.5, 0.3, 1.0);
  b.add(0.5, 0.3, 2.0);
  b.add(0.05, 0.9, 4.0);
  a.merge(b);
  EXPECT_NEAR(std::abs(a.total_weight() - std::complex<double>(7.0, 0.0)), 0.0, 1e-14);
  double at_03 = 0.0, at_09 = 0.0;
  for (const auto &c : a.cells()) {
    if (c.sigma == a.sigma(a.sigma_index(0.3))) { at_03 += c.weight.real(); }
    if (c.sigma == a.sigma(a.sigma_index(0.9))) { at_09 += c.weight.real(); }
  }
  EXPECT_NEAR(at_03, 3.0, 1e-14);
  EXPECT_NEAR(at_09, 4.0, 1e-14);
}

TEST(AdaptiveBins, save_abin) { // NOLINT
  Params P;
  setup(P);
  AdaptiveBinsPosNeg<double> ab(P);
  ab.add(0.5, 0.3, std::complex<double>(2.0, 0.5));
  ab.add(-0.25, 0.7, std::complex<double>(1.0, -0.5));
  const std::string fn = "adaptive_bins_test.abin";
  for (const bool with_imag : {false, true}) {
    ab.save(fn, with_imag);
    std::ifstream F(fn, std::ios::binary);
    std::vector<double> v;
    double x{};
    while (F.read(reinterpret_cast<char *>(&x), sizeof(double))) v.push_back(x);
    const std::size_t width = with_imag ? 4 : 3;
    ASSERT_EQ(v.size() % width, 0);
    const auto records = v.size() / width;
    ASSERT_EQ(records, ab.pos.nr_cells() + ab.neg.nr_cells());
    double re = 0.0, im = 0.0;
    bool seen_negative = false;
    for (std::size_t r = 0; r < records; r++) {
      const auto e = v[r * width], s = v[r * width + 1];
      if (e < 0.0) {
        seen_negative = true;
      } else {
        EXPECT_FALSE(seen_negative) << "positive frequencies come first";
      }
      EXPECT_DOUBLE_EQ(s, e > 0.0 ? 0.3 : 0.7);
      re += v[r * width + 2];
      if (with_imag) im += v[r * width + 3];
    }
    EXPECT_NEAR(re, 3.0, 1e-14);
    if (with_imag) {
      EXPECT_NEAR(im, 0.0, 1e-14);
    }
  }
  std::remove(fn.c_str());
}

TEST(AdaptiveBins, merge_weighted) { // NOLINT
  Params P;
  setup(P);
  AdaptiveBins<double> a(P), b(P);
  b.add(0.5, 0.3, 2.0);
  b.add(0.005, 0.3, 4.0);
  a.merge_weighted(b, [](const double energy) { return energy > 0.1 ? 0.5 : 0.0; });
  EXPECT_NEAR(std::abs(a.total_weight() - std::complex<double>(1.0, 0.0)), 0.0, 1e-14); // only the cells near 0.5, halved
}

TEST(AdaptiveSigma, formula) { // NOLINT
  Params P;
  P.Lambda = 2.0;
  NRG::Eigen<double> ei(std::vector<double>{0.0, 1.0}, 1.0), ej(std::vector<double>{0.5, 2.0}, 1.0);
  {
    const AdaptiveSigma<double> sigma(ei, ej, P); // no dE/dz: the conventional width
    EXPECT_DOUBLE_EQ(sigma(0, 1, 0.3), std::log(2.0));
  }
  ei.values.set_scale(0.25);
  ej.values.set_scale(0.25);
  ei.values.set_dEdz({0.0, -0.7});
  ej.values.set_dEdz({-0.35, -1.4});
  const AdaptiveSigma<double> sigma(ei, ej, P);
  // omega = scale (e_j - e_i) = 0.25 * (2.0 - 1.0); sigma = scale |dEdz_j - dEdz_i| / |omega| = |-1.4 + 0.7| / 1.0
  EXPECT_DOUBLE_EQ(sigma(1, 1, 0.25), 0.7);
  EXPECT_DOUBLE_EQ(sigma(1, 1, -0.25), 0.7);                                   // |omega|
  EXPECT_DOUBLE_EQ(sigma(0, 0, 0.125), 0.25 * 0.35 / 0.125);
  EXPECT_EQ(sigma(0, 0, 0.0), std::numeric_limits<double>::infinity());       // omega = 0: the last sigma bin
}

TEST(AdaptiveBins, rejects_default_limits) { // NOLINT
  Params P; // emin, emax default to -1
  EXPECT_THROW(AdaptiveBins<double>{P}, std::invalid_argument);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS(); // NOLINT
}
