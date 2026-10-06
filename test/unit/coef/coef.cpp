#include <string>
#include <sstream>
#include <stdexcept>
#include <gtest/gtest.h>

#include "test_common.hpp"
#include "test_clean.hpp"
#include <coef.hpp>

using namespace NRG;

TEST(Coef, parse) { // NOLINT
  Params P;
  auto SymSP = setup_Sym<double>(P);
  auto Sym = SymSP.get();
  auto diag = setup_diag_clean<double>(P, Sym);
  Coef<double> coef(P);
  std::string str = 
    "1\n"
    "0.54528747084262258072\n"
    "0.41550946829175445321\n"
    "1\n"
    "0.e-999\n"
    "0.e-998\n";
  std::istringstream ss(str);
  EXPECT_EQ(P.coefchannels, 1);
  coef.xi.read(ss, P.coefchannels);
  coef.zeta.read(ss, P.coefchannels);

  EXPECT_EQ(coef.xi.nr_tabs(), 1);
  EXPECT_EQ(coef.zeta.nr_tabs(), 1);
  EXPECT_EQ(coef.xi.max(0), 1);
  EXPECT_EQ(coef.zeta.max(0), 1);

  EXPECT_EQ(coef.xi(0,0), 0.54528747084262258072);
  EXPECT_EQ(coef.xi(1,0), 0.41550946829175445321); // first index: N, second index: alpha (channel number)
  EXPECT_EQ(coef.zeta(0,0), 0);
  EXPECT_EQ(coef.zeta(1,0), 0);
  
  coef.xi.set(0, 0, 42.);
  EXPECT_EQ(coef.xi(0,0), 42.);
}

TEST(Coef, difference_quotient) { // NOLINT
  // Binary fractions, so that the expected quotients are exact.
  std::istringstream base_ss("1\n1.0\n0.5\n1\n-2.0\n4.0\n");    // two tables, max index 1
  std::istringstream pert_ss("1\n1.125\n0.4375\n1\n-2.0\n3.0\n");
  set_of_coef_tables<double> base, perturbed, derivative;
  base.read(base_ss, 2);
  perturbed.read(pert_ss, 2);
  derivative.set_difference_quotient(base, perturbed, -0.5);  // backward step, z' = z - 0.5
  EXPECT_EQ(derivative.nr_tabs(), 2);
  EXPECT_EQ(derivative.max(0), 1);
  EXPECT_EQ(derivative.max(1), 1);
  EXPECT_EQ(derivative(0, 0), -0.25);
  EXPECT_EQ(derivative(1, 0), 0.125);
  EXPECT_EQ(derivative(0, 1), 0.0);
  EXPECT_EQ(derivative(1, 1), 2.0);
}

TEST(Coef, truncated_matrix_table_header_throws) { // NOLINT
  std::istringstream ss("1 2");
  EXPECT_THROW(read_matrix_table<double>(ss), std::runtime_error);
}

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS(); // NOLINT
}
