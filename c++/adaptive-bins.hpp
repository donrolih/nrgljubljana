// adaptive-bins.hpp - Two-dimensional binning of the spectral data for adaptive broadening
//
// Adaptive broadening (Lee and Weichselbaum, PRB 94, 235127 (2016), Eq. 14) gives each spectral weight its own
// broadening width sigma, so the weights are collected on a grid in both the frequency and sigma.

#ifndef _adaptive_bins_hpp_
#define _adaptive_bins_hpp_

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

#include "traits.hpp"
#include "params.hpp"
#include "basicio.hpp"  // safe_open
#include "numerics.hpp" // reim

namespace NRG {

// Spectral weights binned in |omega| and in the broadening width sigma. The frequency grid is that of Bins::loggrid_std
// (bins.hpp): P.bins points per decade from 10^floor(log10 emin) to 10^ceil(log10 emax), and a weight between two
// points is split between them linearly in log10 |omega|, as in Bins::add_std. The sigma grid is that of the paper:
// sigma_k = (k+1) step for k = 0..nsigma-1, the last point at sigma_max; a weight goes to the nearest point, values
// beyond either end (and non-finite ones) to the end point. Only cells which receive weight are stored.
template<scalar S>
class AdaptiveBins {
 public:
   using t_weight = weight_traits<S>;
   struct Cell {
     double energy;
     double sigma;
     t_weight weight;
   };

 private:
   const Params &P;
   double emin{}, emax{};
   double log10emin{}, log10emax{};
   std::size_t nomega{}; // number of points of the frequency grid
   double sigma_step{};
   std::size_t nsigma{}; // number of points of the sigma grid
   std::unordered_map<std::uint64_t, t_weight> cells_; // key: omega index * nsigma + sigma index

   inline static const double base = 10;

   void add_cell(const std::size_t iomega, const std::size_t isigma, const t_weight weight) {
     cells_[static_cast<std::uint64_t>(iomega) * nsigma + isigma] += weight;
   }

 public:
   explicit AdaptiveBins(const Params &P_) : P(P_) {
     if (P.bins.value() == 0) throw std::invalid_argument("bins must be greater than 0.");
     // Explicit limits are required with adaptive_broadening (Params::validate_adaptive_broadening).
     emin = P.emin;
     emax = P.emax;
     if (!std::isfinite(emin) || !std::isfinite(emax) || emin <= 0.0 || emax <= 0.0 || emin >= emax)
       throw std::invalid_argument("Invalid binning limits: require finite 0 < emin < emax.");
     // As in Bins::setlimits and Bins::loggrid_std.
     log10emin = std::floor(std::log10(emin));
     log10emax = std::ceil(std::log10(emax));
     nomega = static_cast<std::size_t>((log10emax - log10emin) * P.bins + 1.0);
     sigma_step = P.adaptive_broadening_sigma_step;
     const auto sigma_max = P.adaptive_sigma_max();
     if (!(sigma_step > 0.0) || !(sigma_max >= sigma_step))
       throw std::invalid_argument("Invalid sigma grid: require 0 < adaptive_broadening_sigma_step <= adaptive_broadening_sigma_max.");
     nsigma = static_cast<std::size_t>(std::llround(sigma_max / sigma_step));
     my_assert(nsigma >= 1);
   }

   [[nodiscard]] auto omega_points() const noexcept { return nomega; }
   [[nodiscard]] auto sigma_points() const noexcept { return nsigma; }
   [[nodiscard]] double omega(const std::size_t i) const { return std::pow(base, log10emin + double(i) / P.bins); }
   [[nodiscard]] double sigma(const std::size_t k) const { return double(k + 1) * sigma_step; }
   [[nodiscard]] auto nr_cells() const noexcept { return cells_.size(); }

   // Index of the nearest sigma point, clamped to the grid.
   [[nodiscard]] std::size_t sigma_index(const double s) const {
     if (!std::isfinite(s)) return nsigma - 1;
     const auto k = std::llround(s / sigma_step) - 1;
     if (k < 0) return 0;
     return std::min(static_cast<std::size_t>(k), nsigma - 1);
   }

   // Add weight at frequency |omega| = energy >= 0 with broadening width s.
   void add(const double energy, const double s, const t_weight weight) {
     if (std::abs(weight) < P.discard_immediately * energy) return; // as Bins::add
     const auto k = sigma_index(s);
     if (energy < emin) { // as Bins::add_std
       add_cell(0, k, weight);
       return;
     }
     const double x = (std::log10(energy) - log10emin) * P.bins;
     const double int_part = std::floor(x);
     if (int_part < 0) {
       add_cell(0, k, weight);
     } else if (const auto index = std::size_t(int_part); index >= nomega - 1) {
       add_cell(nomega - 1, k, weight);
     } else {
       const double rem = x - int_part;
       add_cell(index, k, (1.0 - rem) * weight);
       add_cell(index + 1, k, rem * weight);
     }
   }

   void merge(const AdaptiveBins<S> &b) {
     my_assert(nomega == b.nomega && nsigma == b.nsigma && log10emin == b.log10emin && sigma_step == b.sigma_step);
     for (const auto &[key, w] : b.cells_) cells_[key] += w;
   }

   void clear() { cells_.clear(); }

   // The stored cells, sorted by frequency and then by sigma.
   [[nodiscard]] std::vector<Cell> cells() const {
     std::vector<std::pair<std::uint64_t, t_weight>> sorted(cells_.begin(), cells_.end());
     std::sort(sorted.begin(), sorted.end(), [](const auto &a, const auto &b) { return a.first < b.first; });
     std::vector<Cell> result;
     result.reserve(sorted.size());
     for (const auto &[key, w] : sorted)
       result.push_back({omega(key / nsigma), sigma(key % nsigma), w});
     return result;
   }

   [[nodiscard]] t_weight total_weight() const {
     t_weight sum{};
     for (const auto &[key, w] : cells_) sum += w;
     return sum;
   }
};

// Adaptive bins for positive and negative frequencies, as ChainBinning (spectrum.hpp) holds Bins.
template<scalar S>
class AdaptiveBinsPosNeg {
 public:
   AdaptiveBins<S> pos, neg; // neg holds |omega|
   explicit AdaptiveBinsPosNeg(const Params &P) : pos(P), neg(P) {}
   void add(const double omega, const double s, const weight_traits<S> weight) {
     if (omega >= 0.0) pos.add(omega, s, weight);
     else neg.add(-omega, s, weight);
   }
   void merge(const AdaptiveBinsPosNeg<S> &b) {
     pos.merge(b.pos);
     neg.merge(b.neg);
   }
   void clear() {
     pos.clear();
     neg.clear();
   }
   // Write the .abin file: headerless native doubles, records (omega, sigma, Re w[, Im w if reim]); positive
   // frequencies first in ascending order, then negative ones from near zero downwards, as savebins writes .bin.
   void save(const std::string &filename, const bool with_imag) const {
     std::ofstream F = safe_open(filename, true); // binary
     const auto write = [&F, with_imag](const double e, const double s, const weight_traits<S> w) {
       const auto [wr, wi] = reim(w);
       F.write(reinterpret_cast<const char *>(&e), sizeof(double));
       F.write(reinterpret_cast<const char *>(&s), sizeof(double));
       F.write(reinterpret_cast<const char *>(&wr), sizeof(double));
       if (with_imag) F.write(reinterpret_cast<const char *>(&wi), sizeof(double));
     };
     for (const auto &c : pos.cells()) write(c.energy, c.sigma, c.weight);
     for (const auto &c : neg.cells()) write(-c.energy, c.sigma, c.weight);
     F.close();
     if (!F) throw std::runtime_error("Error writing " + filename);
   }
};

} // namespace NRG

#endif
