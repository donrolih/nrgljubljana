#include <gtest/gtest.h>

#include <cmath>
#include <functional>
#include <stdexcept>
#include <utility>

#include <mixchain/dz_chain.hpp>
#include <mixchain/star.hpp>

using namespace NRG::MixChain;

namespace {

constexpr double lambda_value = 2.0;
constexpr unsigned int mmax   = 40;
constexpr unsigned int nmax   = 12;
constexpr double z            = 0.5;
constexpr double dz           = 1e-4;

GammaInput<double> input_of(const std::function<Matrix<double>(double)> &gamma) {
  GammaInput<double> input;
  for (int k = 0; k <= 100; k++) {
    input.pos.omega.push_back(0.01 * k);
    input.pos.gamma.push_back(gamma(0.01 * k));
    input.neg.omega.push_back(0.01 * k);
    input.neg.gamma.push_back(gamma(0.01 * k));
  }
  input.channels = static_cast<int>(input.pos.gamma.front().rows());
  return input;
}

Matrix<double> diagonal_of(const double first, const double second) {
  Matrix<double> m = Matrix<double>::Zero(2, 2);
  m(0, 0)          = first;
  m(1, 1)          = second;
  return m;
}

// Two bands of different shape, mixed by an angle that depends on the frequency, as in the chain tests.
Matrix<double> turning(const double omega) {
  const auto angle = 0.4 + 0.8 * omega;
  Matrix<double> u(2, 2);
  u << std::cos(angle), -std::sin(angle), std::sin(angle), std::cos(angle);
  return u * diagonal_of(0.7 - 0.2 * omega, 0.2 + 0.3 * omega * omega) * u.transpose();
}

// The chains at z and z + dz, from one setup of the star stage, as mixchain builds them.
std::pair<Chain<double>, Chain<double>> chains_of(const std::function<Matrix<double>(double)> &gamma,
                                                  const bool split = true) {
  StarOptions options;
  options.Lambda       = NRG::Tools::LambdaCache(lambda_value);
  options.mMAX         = mmax;
  options.split_blocks = split;
  StarDiscretizer<double> discretizer(input_of(gamma), options);
  ChainOptions chain_options;
  chain_options.Nmax = nmax;
  return {build_chain(discretizer.star(z), chain_options), build_chain(discretizer.star(z + dz), chain_options)};
}

} // namespace

TEST(MixChainDzChain, a_turning_gamma_changes_smoothly) { // NOLINT
  const auto [main, perturbed] = chains_of(turning);
  ASSERT_EQ(main.blocks.size(), 1U);
  const auto comparison = compare_dz_chain(main, perturbed, dz);
  // At the fixed point every coefficient scales as Lambda^(-z), a rate of ln(Lambda); the first sites change less.
  const auto log_lambda = std::log(lambda_value);
  EXPECT_GT(comparison.max_hopping_rate, 0.5 * log_lambda);
  EXPECT_LT(comparison.max_hopping_rate, 2.0 * log_lambda);
  EXPECT_LT(comparison.max_onsite_rate, 2.0 * log_lambda);
  // Theta follows the integral of Gamma only approximately where the eigenvectors turn, so it moves with z, but by
  // much less than the chain.
  EXPECT_GT(comparison.theta_change, 0.0);
  EXPECT_LT(comparison.theta_change, 1e-4);
}

TEST(MixChainDzChain, a_diagonal_gamma_keeps_theta) { // NOLINT
  // Each block is a scalar band, whose star reproduces the integral of Gamma exactly: Theta changes only through the
  // weight below the lowest interval, Lambda^(-z-mMAX).
  const auto [main, perturbed] =
    chains_of([](const double omega) { return diagonal_of(0.5, 0.3 + 0.2 * omega * omega); });
  ASSERT_EQ(main.blocks.size(), 2U);
  const auto comparison = compare_dz_chain(main, perturbed, dz);
  EXPECT_LT(comparison.theta_change, 1e-14);
  EXPECT_LT(comparison.max_hopping_rate, 2.0 * std::log(lambda_value));
}

TEST(MixChainDzChain, reports_a_jump) { // NOLINT
  // A hopping that flips its sign between z and z + dz, as a crossing of the gauge would, changes by 2/dz per unit z.
  auto [main, perturbed] = chains_of(turning);
  perturbed.T[5]         = -main.T[5];
  const auto comparison  = compare_dz_chain(main, perturbed, dz);
  EXPECT_EQ(comparison.max_hopping_site, 5U);
  EXPECT_NEAR(comparison.max_hopping_rate, 2.0 / dz, 1e-9 / dz);
}

TEST(MixChainDzChain, a_zero_hopping_has_nothing_to_compare) { // NOLINT
  auto [main, perturbed] = chains_of(turning);
  main.T[nmax].setZero();
  main.E[nmax].setZero();
  perturbed.T[nmax].setZero();
  perturbed.E[nmax].setZero();
  const auto comparison = compare_dz_chain(main, perturbed, dz);
  EXPECT_TRUE(std::isfinite(comparison.max_hopping_rate));
  EXPECT_TRUE(std::isfinite(comparison.max_onsite_rate));
}

TEST(MixChainDzChain, refuses_a_change_of_rank) { // NOLINT
  const auto [main, perturbed] = chains_of(turning);
  auto hopping                 = perturbed;
  hopping.block_diagnostics[0].hopping_ranks[3]--;
  EXPECT_THROW(compare_dz_chain(main, hopping, dz), std::runtime_error);
  auto theta = perturbed;
  theta.block_diagnostics[0].theta_rank--;
  EXPECT_THROW(compare_dz_chain(main, theta, dz), std::runtime_error);
}

TEST(MixChainDzChain, refuses_chains_of_different_structure) { // NOLINT
  const auto [main, perturbed] = chains_of(turning);
  auto shorter                 = perturbed;
  shorter.Nmax--;
  EXPECT_THROW(compare_dz_chain(main, shorter, dz), std::invalid_argument);
  auto other_gauge  = perturbed;
  other_gauge.gauge = ChainGauge::nambu;
  EXPECT_THROW(compare_dz_chain(main, other_gauge, dz), std::invalid_argument);
  EXPECT_THROW(compare_dz_chain(main, perturbed, 0.0), std::invalid_argument);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
