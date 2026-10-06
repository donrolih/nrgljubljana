// broaden_common.hpp - Helpers shared by the broadening tools broaden and abroaden: writing the result, the frequency
// mesh and reading binary spectral files. The kernels and the convolution come from tools/broaden/broadening_core.hpp,
// the parallel loop from parallel.hpp.

#ifndef _common_broaden_common_hpp_
#define _common_broaden_common_hpp_

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <atomic>
#include <chrono>
#include <cstddef>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <cassert>
#include <cfloat>
#include <exception>
#include <utility>
#include <vector>
#include <map>
#include <string>
#include <algorithm>
#include <functional>
#include <ios>
#include <istream>
#include <limits>
#include <mutex>
#include <optional>
#include <ostream>
#include <stdexcept>
#include <numeric>
#include <thread>

#include <range/v3/all.hpp>

#include "../broaden/broadening_core.hpp"
#include "misc.hpp"
#include "basicio.hpp"

#include "diagnostics.hpp"
#include "gsl_config.hpp"
#include "output_file.hpp"
#include "parallel.hpp"

namespace NRG::Broaden {

using std::abs; // important!!
using std::size_t;
using NRG::Tools::for_each_index;
using Core::convolve;
using Core::gaussian_kernel;
using Core::derfd_kernel;

inline std::string tostring(const int i) {
  std::ostringstream S;
  S << i;
  return S.str();
}

template<typename S, typename T>
void save(const std::string &filename, const std::vector<S> &x, const std::vector<T> &y, 
          bool verbose, const int SAVE_PREC = 18) {
  if (verbose) { std::cerr << "Saving " << filename << std::endl; }
  std::ofstream F(filename.c_str());
  if (!F) throw std::runtime_error("Failed to open " + filename + " for writing.");
  F << std::setprecision(SAVE_PREC);
  assert(x.size() == y.size());
  const auto nr = x.size();
  for (auto i = 0; i < nr; i++) F << x[i] << " " << y[i] << '\n';
  NRG::Tools::finish_output(F, filename);
  try {
    F.close();
  } catch (const std::ios_base::failure &error) {
    throw std::runtime_error(filename + ": output close failed: " + error.what());
  }
  if (!F) throw std::runtime_error(filename + ": output close failed.");
}

// Estimate the weight using the trapezoidal rule. The x array must be sorted.
template<typename T> T trapez(const std::vector<T> &x, const std::vector<T> &y) {
  T weight = 0.0;
  assert(x.size() == y.size());
  const auto nr = x.size();
  for (auto i = 1; i < nr; i++) {
    assert(x[i] >= x[i-1]);
    weight += (y[i-1] + y[i]) / 2 * (x[i] - x[i-1]);
  }
  return weight;
}

// Create a mesh on which the output spectral function will be computed. a is the accumulation point of the mesh.
auto make_mesh(const double min, const double max, const double ratio, 
               const double a, const bool add_positive = true, const bool add_negative = false) {
  if (!std::isfinite(ratio) || ratio <= 1.0)
    throw std::invalid_argument("broaden_ratio must be finite and greater than 1.");
  if (!std::isfinite(max) || max <= 0.0)
    throw std::invalid_argument("broaden_max must be finite and greater than 0.");
  if (!std::isfinite(min) || min <= 0.0)
    throw std::invalid_argument("broaden_min must be finite and greater than 0.");
  if (min >= max) throw std::invalid_argument("broaden_min must be smaller than broaden_max.");
  const auto rescale_factor = (max-a)/max;
  std::vector<double> mesh;
  for (double z = max; z > min;) {
    const auto x = a + z * rescale_factor;
    if (add_positive) mesh.push_back(x);
    if (add_negative) mesh.push_back(-x);
    const double next = z / ratio;
    if (!(next < z)) throw std::runtime_error("Frequency mesh generation failed to make progress.");
    z = next;
  }
  std::sort(mesh.begin(), mesh.end());
  return mesh;
}

auto load_mesh(const std::string &filename, bool verbose = false) {
  auto l = readtable<double,double>(filename, false);
  std::vector<double> mesh;
  for (const auto &first : l | ranges::views::transform([](const auto& pair) { return pair.first; }))
    mesh.push_back(first);
  if (verbose) std::cerr << "Reading " << filename << "; " << mesh.size() << " lines read." << std::endl;
  return mesh;
}

auto file_size(std::ifstream &f) {
  f.seekg(0, std::ios::beg);
  const auto begin_pos = f.tellg();
  f.seekg(0, std::ios::end);
  const auto end_pos   = f.tellg();
  return end_pos - begin_pos;
}
   
auto nr_doubles(std::ifstream &f) {
  const auto len = file_size(f);
  if (len % sizeof(double) != 0) throw std::runtime_error("Binary input size is not a multiple of sizeof(double).");
  return len/sizeof(double);
}
  
auto nr_rows(std::ifstream &f, const int cols) {
  const auto d = nr_doubles(f);
  if (d % cols != 0) throw std::runtime_error("Binary input does not contain a whole number of rows.");
  return d / cols;
}
   
// Load a file containing binary representation of raw spectral density. The grid is not assumed to be uniform.
auto load(const std::string &filename, const int nrcol, const bool verbose) {
  std::ifstream f(filename.c_str(), std::ios::in | std::ios::binary);
  if (!f.good() || f.eof() || !f.is_open()) 
    throw std::runtime_error("Error opening file " + filename + " for reading.");
  if (verbose) { std::cerr << "Reading " << filename << std::endl; }
  const auto len = file_size(f); // in bytes
  const auto cols = 1 + nrcol;
  const auto rows = nr_rows(f, cols);
  if (verbose) { std::cerr << "len=" << len << ", " << rows << " data points" << std::endl; }
  auto buffer = std::vector<double>(cols*rows);
  f.seekg(0, std::ios::beg); // Return to the beginning of the file.
  f.read((char *)buffer.data(), len);
  if (f.fail()) throw std::runtime_error("Error reading " + filename);
  f.close();
  return buffer;
}

} // namespace

#endif
