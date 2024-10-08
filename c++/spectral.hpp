// spectral.h - Code for handling spectral function data
// Copyright (C) 2005-2020 Rok Zitko

#ifndef _spectral_hpp_
#define _spectral_hpp_

#include <utility>
#include <vector>
#include <cmath>
#include <limits> // quiet_NaN

#include <range/v3/all.hpp>

#include "traits.hpp"
#include "io.hpp"

namespace NRG {

// Container for holding spectral information represented by delta peaks. "Weight" is of type t_weight (complex).
template<scalar S>
using t_delta_peak = std::pair<double, weight_traits<S>>;

template<scalar S, typename t_weight = weight_traits<S>>
class Spikes : public std::vector<t_delta_peak<S>> {
 public:
   template<typename T>
     void save(T&& F, const int prec, const bool imagpart) {
       F << std::setprecision(prec);
       for (const auto &[e, w] : *this) outputxy(F, e, w, imagpart);
     }
   [[nodiscard]] auto sum_weights() const { return sum2(*this); }
   template<typename F>
     [[nodiscard]] auto sum(const bool invert, F && f) const {
       return ranges::accumulate(*this, t_weight{}, {}, [&f,invert](const auto &x){ const auto &[e,w] = x; return w*f(invert ? -e : e); });
     }
};

#ifndef M_SQRTPI
#define M_SQRTPI 1.7724538509055160273
#endif

// cf. A. Weichselbaum and J. von Delft, cond-mat/0607497

// Modified log-Gaussian broadening kernel. For gamma=alpha/4, the
// kernel is symmetric in both arguments.
inline double BR_L(double e, double ept, double alpha) {
  if ((e < 0.0 && ept > 0.0) || (e > 0.0 && ept < 0.0)) return 0.0;
  if (ept == 0.0) return 0.0;
  const double gamma = alpha/4.0;
  return exp(-pow(log(e / ept) / alpha - gamma, 2)) / (alpha * abs(e) * M_SQRTPI);
}

// Normalized to 1, width omega0. The kernel is symmetric in both
// arguments.
inline double BR_G(double e, double ept, double omega0) { return exp(-pow((e - ept) / omega0, 2)) / (omega0 * M_SQRTPI); }

// Note: 'ept' is the energy of the delta peak in the raw spectrum,
// 'e' is the energy of the data point in the broadened spectrum.
inline double BR_NEW(double e, double ept, double alpha, double omega0) {
  double part_l = BR_L(e, ept, alpha);
  // Most of the time we only need to compute part_l (loggaussian)
  if (abs(e) > omega0) return part_l;
  // Note: this is DIFFERENT from the broadening kernel proposed by
  // by Weichselbaum et al. This BR_h is a function of 'e', not
  // of 'ept'! This breaks the normalization at finite temperatures.
  // On the other hand, it gives nicer spectra when used in conjunction
  // with the self-energy trick.
  double BR_h = exp(-pow(log(abs(e) / omega0) / alpha, 2));
  my_assert(BR_h >= 0.0 && BR_h <= 1.0);
  return part_l * BR_h + BR_G(e, ept, omega0) * (1.0 - BR_h);
}

// Calculate "moment"-th spectral moment.
template<scalar S, typename t_weight = weight_traits<S>>
auto moment(const Spikes<S> &s_neg, const Spikes<S> &s_pos, const int moment) {
  auto sumA = ranges::accumulate(s_pos, t_weight{}, {}, [moment](const auto &x){ const auto &[e,w] = x; return w*pow(e,moment); });
  auto sumB = ranges::accumulate(s_neg, t_weight{}, {}, [moment](const auto &x){ const auto &[e,w] = x; return w*pow(-e,moment); });
  return sumA+sumB;
}

inline constexpr double fermi_fnc(const double omega, const double T) {
  return 1 / (1 + exp(-omega / T));
}

inline constexpr double bose_fnc(const double omega, const double T) {
  const auto d = 1.0 - exp(-omega / T);
  return d != 0.0 ? 1.0/d : std::numeric_limits<double>::quiet_NaN();
}

// Integrated spectral function with a kernel as in FDT for fermions
template<scalar S>
auto fd_fermi(const Spikes<S> &s_neg, const Spikes<S> &s_pos, const double T) {
  const auto fnc = [T](const auto x) { return fermi_fnc(x, T); };
  return s_neg.sum(true, fnc) + s_pos.sum(false, fnc);
}

// Ditto for bosons
template<scalar S>
auto fd_bose(const Spikes<S> &s_neg, const Spikes<S> &s_pos, const double T) {
  const auto fnc = [T](const auto x) { return bose_fnc(x, T); };
  return s_neg.sum(true, fnc) + s_pos.sum(false, fnc);
}

} // namespace

#endif // _spectral_hpp_
