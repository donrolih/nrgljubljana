// Reading the D block of the data file: the Wilson chain at a shifted z, used by adaptive_broadening.

#include <cstdio>
#include <fstream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <gtest/gtest.h>

#include <params.hpp>
#include <read-input.hpp>

using namespace NRG;

namespace {

// data_qs.txt is test/c++/test1_qs/data: QS, one channel, Nmax=4, z block ending the file.
std::string base_data() {
  std::ifstream F("data_qs.txt");
  if (!F) throw std::runtime_error("cannot open data_qs.txt");
  std::ostringstream ss;
  ss << F.rdbuf();
  return ss.str();
}

// A D block with the layout of the z block of data_qs.txt (max index 4, one coefficient channel).
// zeta0 is the on-site energy of site 0 at z', which is part of the initial cluster.
std::string d_block(const std::string &header = "D -0.0025", const size_t maxindex = 4, const double zeta0 = 0.0) {
  std::ostringstream ss;
  ss << header << "\n" << maxindex << "\n";
  for (size_t n = 0; n <= maxindex; n++) ss << 0.5 - 0.1 * double(n) << "\n";
  ss << maxindex << "\n";
  for (size_t n = 0; n <= maxindex; n++) ss << (n == 0 ? zeta0 : 0.0) << "\n";
  return ss.str();
}

std::unique_ptr<Params> make_params(const bool adaptive) {
  const auto filename = "perturbed_chain.param";
  {
    std::ofstream file(filename);
    file << "[param]\n"
         << "symtype=QS\n"
         << "Lambda=2.0\n"
         << "Nmax=4\n"
         << "silent=true\n"
         << "adaptive_broadening=" << (adaptive ? "true" : "false") << "\n";
  }
  auto P = std::make_unique<Params>(filename, "param", std::make_unique<Workdir>(".", true), true, true);
  std::remove(filename);
  return P;
}

// InputData keeps a reference to Params (through Coef), so both are kept together.
struct Run {
  std::unique_ptr<Params> P;
  std::unique_ptr<InputData<double>> input;
};

Run read_data(const bool adaptive, const std::string &extra) {
  Run run{make_params(adaptive), nullptr};
  std::istringstream ss(base_data() + extra);
  run.input = std::make_unique<InputData<double>>(*run.P, ss);
  return run;
}

} // namespace

TEST(PerturbedChain, absent_and_switch_off) { // NOLINT
  const auto run = read_data(false, "");
  EXPECT_FALSE(run.input->dz.has_value());
}

TEST(PerturbedChain, read_with_switch_on) { // NOLINT
  const auto run = read_data(true, d_block());
  const auto &input = *run.input;
  ASSERT_TRUE(input.dz.has_value());
  EXPECT_EQ(input.dz.value(), -0.0025);
  EXPECT_EQ(input.coef_dz.xi.nr_tabs(), 1);
  EXPECT_EQ(input.coef_dz.xi.max(0), 4);
  EXPECT_DOUBLE_EQ(input.coef_dz.xi(1, 0), 0.4);
  EXPECT_DOUBLE_EQ(input.coef_dz.zeta(4, 0), 0.0);
  EXPECT_DOUBLE_EQ(input.coef.xi(0, 0), 0.54528747084262258072); // the z block is untouched
  // The derivative (xi(z') - xi(z))/dz, d_block() has xi(z') = 0.5 - 0.1 n and zeta(z') = 0
  EXPECT_EQ(input.dcoef.xi.nr_tabs(), 1);
  EXPECT_EQ(input.dcoef.xi.max(0), 4);
  EXPECT_DOUBLE_EQ(input.dcoef.xi(0, 0), (0.5 - 0.54528747084262258072) / -0.0025);
  EXPECT_DOUBLE_EQ(input.dcoef.xi(4, 0), (0.1 - 0.17491491876878190148) / -0.0025);
  EXPECT_DOUBLE_EQ(input.dcoef.zeta(2, 0), 0.0);
}

TEST(PerturbedChain, no_derivative_with_switch_off) { // NOLINT
  const auto run = read_data(false, d_block());
  EXPECT_EQ(run.input->dcoef.xi.nr_tabs(), 0);
  EXPECT_EQ(run.input->dcoef.zeta.nr_tabs(), 0);
}

TEST(PerturbedChain, read_but_unused_with_switch_off) { // NOLINT
  const auto run = read_data(false, d_block("D 0.001"));
  const auto &input = *run.input;
  ASSERT_TRUE(input.dz.has_value());
  EXPECT_EQ(input.dz.value(), 0.001);
}

TEST(PerturbedChain, switch_on_requires_block) { // NOLINT
  EXPECT_THROW(read_data(true, ""), std::invalid_argument);
}

TEST(PerturbedChain, rejects_bad_dz) { // NOLINT
  EXPECT_THROW(read_data(true, d_block("D 0")), std::invalid_argument);
  EXPECT_THROW(read_data(true, d_block("D")), std::invalid_argument);
  EXPECT_THROW(read_data(true, d_block("D abc")), std::invalid_argument);
  EXPECT_THROW(read_data(true, d_block("D inf")), std::invalid_argument);
}

TEST(PerturbedChain, rejects_second_block) { // NOLINT
  EXPECT_THROW(read_data(true, d_block() + d_block()), std::invalid_argument);
}

TEST(PerturbedChain, rejects_shape_mismatch) { // NOLINT
  EXPECT_THROW(read_data(true, d_block("D -0.0025", 3)), std::invalid_argument);
}

TEST(PerturbedChain, rejects_z_dependent_initial_cluster) { // NOLINT
  EXPECT_THROW(read_data(true, d_block("D -0.0025", 4, 0.01)), std::invalid_argument);
  EXPECT_NO_THROW(read_data(true, d_block("D -0.0025", 4, 1e-17)));  // round-off level is tolerated
  EXPECT_NO_THROW(read_data(false, d_block("D -0.0025", 4, 0.01))); // irrelevant with the switch off
}

TEST(PerturbedChain, rejects_superconducting_chain) { // NOLINT
  const std::string z_block = "Z\n4\n0\n0\n0\n0\n0\n4\n0\n0\n0\n0\n0\n";
  EXPECT_THROW(read_data(true, d_block() + z_block), std::invalid_argument);
  EXPECT_NO_THROW(read_data(false, d_block() + z_block));
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS(); // NOLINT
}
