// Channel-mixing discretization for NRG
// ** The chain at z and at z' = z + dz, for adaptive broadening

#ifndef _mixchain_dz_chain_hpp_
#define _mixchain_dz_chain_hpp_

#include <algorithm>
#include <cstddef>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <string>

#include "blocks.hpp"
#include "chain.hpp"
#include "types.hpp"

namespace NRG::MixChain {

// How the chain at z' = z + dz differs from the chain at z. nrg forms the derivative of the chain in z from the
// difference of the two, which means something only where the chain is a smooth function of z between them: in the
// same gauge (both are built with the same options), with the same rank structure, and with differences of order dz.
// The rates are relative changes per unit z; at the fixed point of a Wilson chain they are of order ln(Lambda).
//
// The impurity coupling V enters the initial cluster, whose derivative nrg does not take. Theta = V V^dagger is the
// sum of the couplings of the star, which reproduces the integral of Gamma exactly in the trace but only approximately
// off the diagonal when the eigenvectors of Gamma turn, so for a matrix Gamma it can depend on z; that change is
// neglected in dH/dz, and theta_change says how large it is.
struct DzComparison {
  double theta_change{};             // ||Theta(z') - Theta(z)|| / ||Theta(z)||
  double max_hopping_rate{};         // max_n ||T_n(z') - T_n(z)|| / (||T_n(z)|| |dz|)
  unsigned int max_hopping_site{};
  double max_onsite_rate{};          // max_n ||E_n(z') - E_n(z)|| / (max(||E_n(z)||, ||T_n(z)||) |dz|)
  unsigned int max_onsite_site{};
};

namespace detail {

inline std::string rank_mismatch(const std::string &what, const int at_z, const int at_perturbed) {
  std::ostringstream message;
  message << "The chain at z+dz has " << what << " " << at_perturbed << ", the chain at z " << at_z
          << ": the derivative of the chain in z is not defined there. A smaller Nmax, or a dz that does not cross the "
          << "change of rank, avoids it.";
  return message.str();
}

} // namespace detail

// The two chains must have the same structure: anything else means they were not built from the same Gamma with the
// same options, or that the rank of the chain changes between z and z', where no derivative exists.
template<typename S>
DzComparison compare_dz_chain(const Chain<S> &main, const Chain<S> &perturbed, const double dz) {
  if (!(std::isfinite(dz) && dz != 0.0)) throw std::invalid_argument("dz must be a nonzero finite number.");
  if (main.channels != perturbed.channels || main.Nmax != perturbed.Nmax || main.gauge != perturbed.gauge
      || main.blocks != perturbed.blocks || main.E.size() != perturbed.E.size()
      || main.T.size() != perturbed.T.size())
    throw std::invalid_argument("The chains at z and z+dz differ in channels, Nmax, gauge or blocks.");
  // Block by block: the ranks of the whole chain are sums over the blocks, which could hide a change in two of them.
  if (main.block_diagnostics.size() != perturbed.block_diagnostics.size())
    throw std::invalid_argument("The chains at z and z+dz differ in their blocks.");
  for (std::size_t b = 0; b < main.block_diagnostics.size(); b++) {
    const auto &at_z         = main.block_diagnostics[b];
    const auto &at_perturbed = perturbed.block_diagnostics[b];
    const auto block         = main.blocks.size() > 1 ? " in block " + blocks_name({main.blocks[b]}) : std::string();
    if (at_z.theta_rank != at_perturbed.theta_rank)
      throw std::runtime_error(detail::rank_mismatch("Theta of rank" + block, at_z.theta_rank, at_perturbed.theta_rank));
    if (at_z.hopping_ranks.size() != at_perturbed.hopping_ranks.size())
      throw std::invalid_argument("The chains at z and z+dz differ in length.");
    for (std::size_t n = 0; n < at_z.hopping_ranks.size(); n++)
      if (at_z.hopping_ranks[n] != at_perturbed.hopping_ranks[n])
        throw std::runtime_error(detail::rank_mismatch("a hopping at site " + std::to_string(n) + block + " of rank",
                                                       at_z.hopping_ranks[n], at_perturbed.hopping_ranks[n]));
  }

  DzComparison result;
  const Matrix<S> theta           = main.V * main.V.adjoint();
  const Matrix<S> theta_perturbed = perturbed.V * perturbed.V.adjoint();
  const auto theta_norm           = theta.norm();
  if (theta_norm > 0.0) result.theta_change = (theta_perturbed - theta).norm() / theta_norm;

  const auto step = std::abs(dz);
  for (std::size_t n = 0; n < main.T.size(); n++) {
    const auto site    = static_cast<unsigned int>(n);
    const auto hopping = main.T[n].norm();
    // A hopping that is zero in both, past a rank drop, has nothing to compare.
    if (hopping > 0.0) {
      const auto rate = (perturbed.T[n] - main.T[n]).norm() / (hopping * step);
      if (rate > result.max_hopping_rate) {
        result.max_hopping_rate = rate;
        result.max_hopping_site = site;
      }
    }
    // The on-site blocks vanish for a particle-hole symmetric Gamma, so they are measured on the scale of the site.
    const auto scale = std::max(main.E[n].norm(), hopping);
    if (scale > 0.0) {
      const auto rate = (perturbed.E[n] - main.E[n]).norm() / (scale * step);
      if (rate > result.max_onsite_rate) {
        result.max_onsite_rate = rate;
        result.max_onsite_site = site;
      }
    }
  }
  return result;
}

} // namespace NRG::MixChain

#endif
