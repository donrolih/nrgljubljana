#ifndef _read_input_hpp_
#define _read_input_hpp_

#include <algorithm>
#include <cctype>
#include <cstddef>
#include <fstream>
#include <ios>
#include <istream>
#include <cmath>
#include <memory>
#include <optional>
#include <ostream>
#include <string>
#include <iostream>
#include <stdexcept>
#include <tuple>
#include <utility>

#include "traits.hpp"
#include "constants.hpp"
#include "symmetry.hpp"
#include "params.hpp"
#include "eigen.hpp"
#include "operators.hpp"
#include "floquet.hpp"
#include "coef.hpp"
#include "tridiag.hpp"
#include "mk_sym.hpp"

#include <fmt/format.h>

namespace NRG {

struct DataHeader {
  std::string symtype;
  size_t channels;
  size_t Nmax;
  size_t nsubs;
};

inline auto read_data_header(std::istream &fdata, const int expected_version = 9) {
  std::string sym_string;
  int dataversion = -1;
  while (fdata.peek() == '#') {
    fdata.ignore();
    if (fdata.peek() == '!') {
      fdata.ignore();
      fdata >> dataversion;
    } else {
      std::string line;
      std::getline(fdata, line);
      if (const auto pos = line.find("symtype", 1); pos != std::string::npos) {
        if (const auto p = line.find_last_of(" \t"); p != std::string::npos && p < line.size() - 1)
          sym_string = line.substr(p + 1);
      }
    }
    if (fdata.peek() == '\n') fdata.ignore();
  }
  my_assert(dataversion == expected_version);
  return DataHeader{sym_string, read_one<size_t>(fdata), read_one<size_t>(fdata), read_one<size_t>(fdata)};
}

// Determine Nmax & Nlen, taking into account P.substeps. Must be called after the
// tridiagonalization routines if not using the tables in the data file.
template<scalar S>
void determine_Nmax_Nlen(const Coef<S> &coef, const size_t Nmax0, Params &P) { // Params is non-const !
  const auto length_coef_table = coef.xi.max(0); // all channels have the same nr. of coefficients
  my_assert(length_coef_table == Nmax0); // check consistency
  P.Nmax = !P.substeps ? Nmax0 : P.channels * Nmax0;
  P.Nlen = P.Nmax;
  if (!P.silent) std::cout << "\nlength_coef_table=" << length_coef_table << " Nmax0=" << Nmax0 << " Nmax=" << P.Nmax << "\n\n";
  my_assert(P.Nlen < MAX_NDX);
}

template<scalar S>
class InputData {
private:
  static auto read_next_block(std::istream &fdata) {
    while (!fdata.eof() && std::isspace(fdata.peek())) fdata.get();
    char ch = char(fdata.get());
    std::string opname;
    std::getline(fdata, opname);
    return std::make_pair(ch, opname);
  }

  void initialize_symmetry(Params &P, const DataHeader &header) {
    my_assert(header.symtype == P.symtype.value());
    Sym = set_symmetry<S>(P, header.symtype, header.channels);
  }

  void read_seed_data(std::istream &fdata, const DataHeader &header, Params &P) {
    diag = DiagInfo<S>(fdata, header.nsubs, P);
    if (!P.silent)
      diag.states_report(Sym->multfnc());
    operators.opch = Opch<S>(fdata, diag, P);
  }

  // D block: the header line holds the signed shift dz = z' - z, followed by the xi and zeta tables of the chain at
  // z', with the layout and units of the z block.
  void read_perturbed_chain(std::istream &fdata, const std::string &header, const Params &P) {
    if (dz) throw std::invalid_argument("The data file contains more than one D block.");
    double value{};
    try {
      value = from_string<double>(strip_whitespace(header));
    } catch (const std::exception &) {
      throw std::invalid_argument(fmt::format("D block: dz must be a finite nonzero number, got [{}].", header));
    }
    if (!std::isfinite(value) || value == 0.0)
      throw std::invalid_argument(fmt::format("D block: dz must be a finite nonzero number, got [{}].", header));
    dz = value;
    coef_dz.xi.read(fdata, P.coefchannels);
    coef_dz.zeta.read(fdata, P.coefchannels);
  }

  // The perturbed chain must have the shape of the chain itself, and only the xi/zeta chains have a perturbed
  // counterpart so far.
  void validate_perturbed_chain(const Params &P) const {
    if (dz && !P.adaptive_broadening && !P.silent)
      std::cout << "The data file contains a D block, which is unused without adaptive_broadening." << std::endl;
    if (!P.adaptive_broadening) return;
    if (!dz) throw std::invalid_argument("adaptive_broadening requires a D block (the chain at a shifted z) in the data file.");
    if (std::string(P.tri) == "cpp")
      throw std::invalid_argument("adaptive_broadening is not supported with tri=cpp.");
    if (coef.delta.nr_tabs() || coef.kappa.nr_tabs())
      throw std::invalid_argument("adaptive_broadening does not support superconducting chains (Z block) yet.");
    if (coef.xiR.nr_tabs() || coef.zetaR.nr_tabs())
      throw std::invalid_argument("adaptive_broadening does not support channel-mixing rungs (X block) yet.");
    const auto check_shape = [](const auto &base, const auto &perturbed, const std::string &name) {
      if (perturbed.nr_tabs() != base.nr_tabs())
        throw std::invalid_argument(fmt::format("D block: {} has {} tables, the z block {}.", name, perturbed.nr_tabs(), base.nr_tabs()));
      for (size_t ch = 0; ch < base.nr_tabs(); ch++)
        if (perturbed.max(ch) != base.max(ch))
          throw std::invalid_argument(fmt::format("D block: {} table {} ends at index {}, the z block at {}.", name, ch, perturbed.max(ch), base.max(ch)));
    };
    check_shape(coef.xi, coef_dz.xi, "xi");
    check_shape(coef.zeta, coef_dz.zeta, "zeta");
    if (P.floquet) throw std::invalid_argument("adaptive_broadening is not supported with floquet.");
    if (P.substeps) throw std::invalid_argument("adaptive_broadening is not supported with substeps.");
    if (operators.ops.count("dHdz"))
      throw std::invalid_argument("adaptive_broadening does not read dH/dz of the initial cluster (singlet dHdz) yet.");
    // The coefficients inside the initial cluster (zeta_0..Ninit, xi_0..Ninit-1) enter the Hamiltonian from the data
    // file, whose derivative is not available yet, so they must not depend on z. Differences at the round-off level
    // of the chain are tolerated.
    double scale = 0.0;
    for (size_t ch = 0; ch < coef.xi.nr_tabs(); ch++)
      for (size_t n = 0; n <= coef.xi.max(ch); n++) scale = std::max(scale, std::abs(coef.xi(n, ch)));
    const auto tolerance = 1e-12 * scale;
    const auto Ninit = size_t(P.Ninit);
    const auto check_initial = [Ninit, tolerance](const auto &base, const auto &perturbed, const std::string &name, const size_t last) {
      for (size_t ch = 0; ch < base.nr_tabs(); ch++)
        for (size_t n = 0; n <= std::min(last, base.max(ch)); n++)
          if (std::abs(perturbed(n, ch) - base(n, ch)) > tolerance)
            throw std::invalid_argument(fmt::format("adaptive_broadening: {}_{} (table {}) depends on z, but it is part of the "
                                                    "initial cluster (Ninit={}), whose dH/dz is not supported yet.", name, n, ch, Ninit));
    };
    check_initial(coef.zeta, coef_dz.zeta, "zeta", Ninit);
    if (Ninit > 0) check_initial(coef.xi, coef_dz.xi, "xi", Ninit - 1);
  }

  void read_blocks(std::istream &fdata, Params &P) {
    while (true) {
      const auto [ch, opname] = read_next_block(fdata);
      if (fdata.eof()) break;
      if (ch != '#' && !P.silent) std::cout << "Reading <||" << opname << "||> (" << ch << ")" << std::endl;
      switch (ch) {
      case '#':
        break;
      case 'e': GS_energy = read_one<double>(fdata); break;
      case 's': operators.ops[opname]  = MatrixElements<S>(fdata, diag); break;
      case 'p': operators.opsp[opname] = MatrixElements<S>(fdata, diag); break;
      case 'g': operators.opsg[opname] = MatrixElements<S>(fdata, diag); break;
      case 'd': operators.opd[opname]  = MatrixElements<S>(fdata, diag); break;
      case 't': operators.opt[opname]  = MatrixElements<S>(fdata, diag); break;
      case 'o': operators.opot[opname] = MatrixElements<S>(fdata, diag); break;
      case 'q': operators.opq[opname]  = MatrixElements<S>(fdata, diag); break;
      case 'z':
        coef.xi.read(fdata, P.coefchannels);
        coef.zeta.read(fdata, P.coefchannels);
        break;
      case 'Z':
        coef.delta.read(fdata, P.coefchannels);
        coef.kappa.read(fdata, P.coefchannels);
        break;
      case 'X':
        coef.xiR.read(fdata, P.coefchannels);
        coef.zetaR.read(fdata, P.coefchannels);
        break;
      case 'D': // Wilson chain at the shifted z' = z + dz, for adaptive broadening; the header line gives dz
        read_perturbed_chain(fdata, opname, P);
        break;
      case 'T':
        coef.ep.read(fdata, P.coefchannels);
        coef.em.read(fdata, P.coefchannels);
        coef.u0p.read(fdata, P.coefchannels);
        coef.u0m.read(fdata, P.coefchannels);
        break;
#ifdef CHAIN_COEF
      case 'W':
        coef.eps.read(fdata);
        coef.t.read(fdata);
        break;
#endif
      default: throw std::invalid_argument(fmt::format("Unknown block {} in data file.", ch));
      }
    }
  }

  // d/dz of the chain coefficients, from the chain at z (z block) and at z + dz (D block)
  void form_chain_derivative() {
    dcoef.xi.set_difference_quotient(coef.xi, coef_dz.xi, dz.value());
    dcoef.zeta.set_difference_quotient(coef.zeta, coef_dz.zeta, dz.value());
  }

  void finalize_coefficients(Params &P, const DataHeader &header) {
    validate_perturbed_chain(P);
    if (P.adaptive_broadening) form_chain_derivative();
    if (std::string(P.tri) == "cpp") Tridiag<S>(coef, header.Nmax, P);
    determine_Nmax_Nlen(coef, header.Nmax, P);
    P.validate_after_data_file();
  }

  void read(std::istream &fdata, Params &P) {
    const auto header = read_data_header(fdata);
    initialize_symmetry(P, header);
    read_seed_data(fdata, header, P);
    read_blocks(fdata, P);
    if (P.floquet)
      P.set_floquet_mode_bounds(validate_floquet_mode_operator(diag, operators));
    finalize_coefficients(P, header);
  }
public:
  std::shared_ptr<Symmetry<S>> Sym;
  DiagInfo<S> diag;
  Operators<S> operators;
  Coef<S> coef;
  Coef<S> coef_dz;          // chain at z' = z + dz (D block), for adaptive broadening
  std::optional<double> dz; // set only if the data file has a D block
  Coef<S> dcoef;            // d/dz of xi and zeta, formed only with adaptive_broadening
  double GS_energy = 0.0;
  InputData(Params &P, std::istream &fdata) : coef(P), coef_dz(P), dcoef(P) {
    read(fdata, P);
  }
  InputData(Params &P, const std::string &filename = "data") : coef(P), coef_dz(P), dcoef(P) {
    std::ifstream fdata(filename);
    if (!fdata) throw std::runtime_error(fmt::format("Can't load initial data from {}.", filename));
    read(fdata, P);
  };
};

} // namespace

#endif
