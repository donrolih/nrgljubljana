#ifndef _truncation_h_
#define _truncation_h_

// Determine the number of states to be retained. Returns Emax - the highest energy to still be retained.
template<typename S>
auto highest_retained_energy(const Step &step, const DiagInfo<S> &diag, const Params &P) {
  auto energies = diag.sorted_energies();
  my_assert(energies.front() == 0.0); // check for the subtraction of Egs
  const auto totalnumber = energies.size();
  size_t nrkeep;
  if (P.keepenergy <= 0.0) {
    nrkeep = P.keep;
  } else {
    const auto keepenergy = P.keepenergy * step.unscale();
    // We add 1 for historical reasons. We thus keep states with E<=Emax, and one additional state which has E>Emax.
    nrkeep = 1 + ranges::count_if(energies, [=](double e) { return e <= keepenergy; });
    nrkeep = std::clamp<size_t>(nrkeep, P.keepmin, P.keep);
  }
  // Check for near degeneracy and ensure that the truncation occurs in a "gap" between nearly-degenerate clusters of
  // eigenvalues.
  if (P.safeguard > 0.0) {
    size_t cnt_extra = 0;
    while (nrkeep < totalnumber && (energies[nrkeep] - energies[nrkeep - 1]) <= P.safeguard && cnt_extra < P.safeguardmax) {
      nrkeep++;
      cnt_extra++;
    }
    if (cnt_extra) debug("Safeguard: keep additional " << cnt_extra << " states");
  }
  nrkeep = std::clamp<size_t>(nrkeep, 1, totalnumber);
  return energies[nrkeep - 1];
}

struct truncate_stats {
  size_t nrall, nrallmult, nrkept, nrkeptmult;
  template<typename S> truncate_stats(const DiagInfo<S> &diag, std::shared_ptr<Symmetry<S>> Sym) {
    nrall = ranges::accumulate(diag, 0, [](int n, const auto &d) { const auto &[I, eig] = d; return n+eig.getdim(); });
    nrallmult = ranges::accumulate(diag, 0, [Sym](int n, const auto &d) { const auto &[I, eig] = d; return n+Sym->mult(I)*eig.getdim(); });
    nrkept = ranges::accumulate(diag, 0, [](int n, const auto &d) { const auto &[I, eig] = d; return n+eig.getnrkept(); });
    nrkeptmult = ranges::accumulate(diag, 0, [Sym](int n, const auto &d) { const auto &[I, eig] = d; return n+Sym->mult(I)*eig.getnrkept(); });
  }
  void report() {
    std::cout << nrgdump4(nrkept, nrkeptmult, nrall, nrallmult) << std::endl;
  }
};

struct NotEnough : public std::exception {};

// Compute the number of states to keep in each subspace. Returns true if an insufficient number of states has been
// obtained in the diagonalization and we need to compute more states.
template<typename S>
void truncate_prepare(const Step &step, DiagInfo<S> &diag, std::shared_ptr<Symmetry<S>> Sym, const Params &P) {
  const auto Emax = highest_retained_energy(step, diag, P);
  for (auto &[I, eig] : diag)
    diag[I].truncate_prepare(step.last() && P.keep_all_states_in_last_step() ? eig.getnrcomputed() :
                             ranges::count_if(eig.value_zero, [Emax](const double e) { return e <= Emax; }));
  std::cout << "Emax=" << Emax/step.unscale() << " ";
  truncate_stats ts(diag, Sym);
  ts.report();
  if (ranges::any_of(diag, [Emax](const auto &d) { const auto &[I, eig] = d; 
    return eig.getnrkept() == eig.getnrcomputed() && eig.value_zero(eig.getnrcomputed()-1) != Emax && eig.getnrcomputed() < eig.getdim(); }))
      throw NotEnough();
  const double ratio = double(ts.nrkept) / ts.nrall;
  fmt::print(FMT_STRING("Kept: {} out of {}, ratio={:.3}\n"), ts.nrkept, ts.nrall, ratio);
}

#endif