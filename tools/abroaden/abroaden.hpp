// abroaden - Adaptive broadening of the binned spectral data of NRG Ljubljana
//
// S.-S. B. Lee and A. Weichselbaum, Phys. Rev. B 94, 235127 (2016). nrg with adaptive_broadening writes .abin files:
// spectral weights binned in the frequency omega and in their broadening width sigma. abroaden averages them over
// the Nz runs (Eq. 15), broadens each cell with a log-Gaussian of width alpha sigma/ln(Lambda) (Eq. 17), and then
// convolves the result over all frequencies with -f'(omega) at the effective temperature gamma (Eq. 16, 21).

#ifndef _abroaden_abroaden_hpp_
#define _abroaden_abroaden_hpp_

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <map>
#include <optional>
#include <ostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <getopt.h>

#include "broadening.hpp"

#include "../common/broaden_common.hpp"
#include "../common/gsl_config.hpp"
#include "../common/parallel.hpp"

namespace NRG::Broaden {

class AdaptiveBroaden {
 private:
   bool verbose = false;
   double broaden_min = 1e-7; // the frequency mesh, as in broaden
   double broaden_max = 2.0;
   double broaden_ratio = 1.01;
   std::string mesh_filename;
   bool one = false;       // Nz=1: read <name> rather than 1/<name>
   bool with_imag = false; // records (omega, sigma, Re w, Im w) rather than (omega, sigma, Re w)
   double gamma_ratio = 0.1; // gamma = gamma_ratio * T (Eq. 21); the paper finds 0.1 typical
   double sigma_min = 0.0;   // optional floor on sigma
   size_t jobs = 1;
   std::string name;
   int Nz = 0;
   double alpha = 0.0;
   double T = 0.0;
   double Lambda = 0.0;
   const std::string output_filename = "spec.dat";

   // A cell of the z-averaged data: frequency, kernel width alpha sigma/ln(Lambda), weight.
   struct Cell {
     double logabs; // ln |omega|, the sort key
     double omega;
     double width;
     double weight;
   };
   std::vector<Cell> positive, negative; // sorted by ln |omega|
   double max_width = 0.0;
   double total_weight = 0.0;
   std::vector<double> mesh;
   std::vector<double> a;

   void usage(std::ostream &F = std::cout) const {
     F << "Usage: abroaden [options] <name> <Nz> <alpha> <T> <Lambda>\n"
       << "\n"
       << "Adaptive broadening (Lee and Weichselbaum, PRB 94, 235127 (2016)) of the .abin files written by nrg with\n"
       << "adaptive_broadening: i/<name> for i=1..Nz are averaged cell by cell, each cell is broadened with a\n"
       << "log-Gaussian of width alpha*sigma/ln(Lambda), and the result is convolved with -f'(omega) at gamma.\n"
       << "alpha is the width of a weight with sigma = ln(Lambda), as in broaden; the paper's alpha corresponds to\n"
       << "alpha*Nz/ln(Lambda).\n"
       << "\n"
       << "Options:\n"
       << " -h -- show help (as sole switch)\n"
       << " -v -- verbose\n"
       << " -V, --version -- show version and exit\n"
       << " -j <n>, --jobs <n> -- worker count (default: first OMP_NUM_THREADS value, or 1)\n"
       << " -m <min>, -M <max>, -r <ratio> -- frequency mesh, as in broaden (default 1e-7, 2, 1.01)\n"
       << " -L <filename> -- load the frequency mesh from a file\n"
       << " -o -- read <name> directly (Nz=1)\n"
       << " -I -- records with an imaginary part (reim=true in nrg); the real part is broadened\n"
       << " -x <ratio> -- gamma = ratio*T for the final -f' convolution (default 0.1)\n"
       << " -S <sigma_min> -- floor on sigma (default none)\n";
   }

   void cmd_line(int argc, char *argv[]) {
     if (argc == 2 && std::string(argv[1]) == "-h") {
       usage();
       std::exit(EXIT_SUCCESS);
     }
     static const option long_options[] = {{"jobs", required_argument, nullptr, 'j'}, {nullptr, 0, nullptr, 0}};
     std::optional<size_t> requested_jobs;
     int c_;
     while ((c_ = getopt_long(argc, argv, "vj:m:M:r:L:oIx:S:", long_options, nullptr)) != -1) {
       switch (c_) {
       case 'v': verbose = true; break;
       case 'j': requested_jobs = NRG::Tools::parse_positive_size(optarg, "Worker count"); break;
       case 'm': broaden_min = std::atof(optarg); break;
       case 'M': broaden_max = std::atof(optarg); break;
       case 'r': broaden_ratio = std::atof(optarg); break;
       case 'L': mesh_filename = optarg; break;
       case 'o': one = true; break;
       case 'I': with_imag = true; break;
       case 'x': gamma_ratio = std::atof(optarg); break;
       case 'S': sigma_min = std::atof(optarg); break;
       default: throw std::invalid_argument("Unknown command-line option.");
       }
     }
     if (requested_jobs) {
       jobs = *requested_jobs;
     } else if (const auto *environment = std::getenv("OMP_NUM_THREADS")) {
       jobs = NRG::Tools::parse_worker_count(environment, "OMP_NUM_THREADS");
     }
     if (argc - optind != 5) {
       usage(std::cerr);
       std::exit(1);
     }
     name = argv[optind];
     Nz = std::atoi(argv[optind + 1]);
     alpha = std::atof(argv[optind + 2]);
     T = std::atof(argv[optind + 3]);
     Lambda = std::atof(argv[optind + 4]);
     if (!(Nz >= 1)) throw std::invalid_argument("Nz must be greater than or equal to 1.");
     if (one && Nz != 1) throw std::invalid_argument("-o requires Nz=1.");
     if (!(alpha > 0.0) || !std::isfinite(alpha)) throw std::invalid_argument("alpha must be finite and greater than 0.");
     if (!(T > 0.0) || !std::isfinite(T)) throw std::invalid_argument("T must be finite and greater than 0.");
     if (!(Lambda > 1.0) || !std::isfinite(Lambda)) throw std::invalid_argument("Lambda must be finite and greater than 1.");
     if (!(gamma_ratio > 0.0) || !std::isfinite(gamma_ratio)) throw std::invalid_argument("The -x ratio must be finite and greater than 0.");
     if (!(sigma_min >= 0.0) || !std::isfinite(sigma_min)) throw std::invalid_argument("sigma_min must be finite and non-negative.");
     std::cout << "Processing: " << name << std::endl;
   }

   // Read i/<name>, i=1..Nz, and average the weights cell by cell (Eq. 15).
   void read_and_merge() {
     const int nrcol = with_imag ? 3 : 2; // columns after omega: sigma, Re w[, Im w]
     const auto cols = 1 + nrcol;
     std::map<std::pair<double, double>, double> cells;
     for (int i = 1; i <= Nz; i++) {
       const auto filename = one ? name : std::to_string(i) + "/" + name;
       const auto buffer = load(filename, nrcol, verbose);
       for (size_t r = 0; r + cols <= buffer.size(); r += cols) {
         const auto omega = buffer[r], sigma = buffer[r + 1], weight = buffer[r + 2];
         if (!std::isfinite(omega) || !std::isfinite(sigma) || !std::isfinite(weight))
           throw std::runtime_error(filename + ": non-finite record.");
         cells[{omega, sigma}] += weight;
       }
     }
     const auto lnL = std::log(Lambda);
     for (const auto &[key, w] : cells) {
       const auto &[omega, sigma] = key;
       const auto weight = w / Nz;
       total_weight += weight;
       if (omega == 0.0 || weight == 0.0) continue;
       const auto width = alpha * std::max(sigma, sigma_min) / lnL;
       if (!(width > 0.0)) throw std::runtime_error("A cell has a non-positive broadening width; use -S.");
       max_width = std::max(max_width, width);
       (omega > 0.0 ? positive : negative).push_back({std::log(std::abs(omega)), omega, width, weight});
     }
     const auto by_log = [](const Cell &x, const Cell &y) { return x.logabs < y.logabs; };
     std::sort(positive.begin(), positive.end(), by_log);
     std::sort(negative.begin(), negative.end(), by_log);
     std::cout << cells.size() << " cells, total weight=" << total_weight << std::endl;
   }

   // Each cell broadened with its own log-Gaussian (Eq. 17a). Cells further than b (7 + b/4) in ln|omega| from the
   // output point, b the largest width, contribute less than exp(-49) relative and are skipped.
   std::vector<double> broaden_log() const {
     std::vector<double> result(mesh.size());
     const auto window = max_width * (7.0 + max_width / 4.0);
     for_each_index(mesh.size(), jobs, [&](const size_t index) {
       const auto m = mesh[index];
       const auto &cells = m > 0.0 ? positive : negative;
       const auto logm = std::log(std::abs(m));
       const auto lo = std::lower_bound(cells.begin(), cells.end(), logm - window,
                                        [](const Cell &c, const double v) { return c.logabs < v; });
       double sum = 0.0;
       for (auto it = lo; it != cells.end() && it->logabs <= logm + window; ++it)
         sum += it->weight * NRG::Broadening::log_gaussian(m, it->omega, it->width);
       result[index] = sum;
     });
     return result;
   }

 public:
   AdaptiveBroaden(int argc, char *argv[]) {
     std::cout << "abroaden - adaptive broadening tool" << std::endl;
     cmd_line(argc, argv);
   }
   void calc() {
     read_and_merge();
     mesh = mesh_filename.empty() ? make_mesh(broaden_min, broaden_max, broaden_ratio, 0.0, true, true)
                                  : load_mesh(mesh_filename, verbose);
     if (mesh.empty()) throw std::runtime_error("Output mesh is empty.");
     if (verbose)
       std::cerr << "alpha=" << alpha << " T=" << T << " Lambda=" << Lambda << " gamma=" << gamma_ratio * T
                 << " sigma_min=" << sigma_min << " mesh points=" << mesh.size() << " jobs=" << jobs << std::endl;
     a = broaden_log();
     convolve(mesh, a, gamma_ratio * T, jobs, derfd_kernel); // Eq. 16: after the log-Gaussian, over all frequencies
     std::cout << "Estimated weight (trapezoidal rule)=" << trapez(mesh, a) << std::endl;
     save(output_filename, mesh, a, verbose);
   }
};

} // namespace NRG::Broaden

#endif
