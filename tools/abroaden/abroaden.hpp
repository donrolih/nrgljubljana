// abroaden - Adaptive broadening of the binned spectral data of NRG Ljubljana
//
// S.-S. B. Lee and A. Weichselbaum, Phys. Rev. B 94, 235127 (2016). nrg with adaptive_broadening writes .abin files:
// spectral weights binned in the frequency omega and in their broadening width sigma. abroaden averages them over
// the Nz runs (Eq. 15), broadens each cell with a log-Gaussian of width alpha sigma/ln(Lambda) (Eq. 17), and then
// convolves the result over all frequencies with -f'(omega) at the effective temperature gamma (Eq. 16, 21).

#ifndef _abroaden_abroaden_hpp_
#define _abroaden_abroaden_hpp_

#include <algorithm>
#include <array>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <map>
#include <optional>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
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
   bool auto_gamma = false;      // choose gamma by the eps_rms scan (Sec. V.2)
   double eps_threshold = 0.002; // threshold on eps_rms for the automatic gamma
   const std::string scan_filename = "gamma_scan.dat";
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
   // The cells of one sign, grouped by kernel width: a group holds the widths in [2^k, 2^(k+1)), sorted by ln|omega|,
   // and the window in ln|omega| beyond which its widest kernel is negligible. A few very wide kernels (large sigma, at
   // |omega| below T) then do not widen the window of all the narrow ones.
   struct Group {
     double window;
     std::vector<Cell> cells;
   };
   std::vector<Group> positive, negative;
   double total_weight = 0.0;
   std::vector<double> mesh;
   std::vector<double> a;
   size_t nr_cells = 0;
   std::chrono::steady_clock::time_point wall_start;
   std::clock_t cpu_start;
   std::vector<std::pair<std::string, double>> phases; // name, wall seconds

   // Run fnc and record its wall time under name.
   template <typename F> auto timed(const std::string &name, F &&fnc) {
     const auto start = std::chrono::steady_clock::now();
     if constexpr (std::is_void_v<decltype(fnc())>) {
       fnc();
       phases.emplace_back(name, std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count());
     } else {
       auto result = fnc();
       phases.emplace_back(name, std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count());
       return result;
     }
   }

   // As broaden's report_timing, plus the phases and the numbers of cells and groups.
   void report_timing() const {
     const auto wall_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - wall_start).count();
     const auto cpu_end = std::clock();
     std::optional<double> cpu_seconds;
     if (cpu_start != static_cast<std::clock_t>(-1) && cpu_end != static_cast<std::clock_t>(-1) && cpu_end >= cpu_start)
       cpu_seconds = static_cast<double>(cpu_end - cpu_start) / static_cast<double>(CLOCKS_PER_SEC);
     std::ostringstream report;
     report << std::setprecision(6) << "Time elapsed: " << wall_seconds << " s\n";
     if (verbose) {
       const auto throughput = wall_seconds > 0.0 ? static_cast<double>(mesh.size()) / wall_seconds : 0.0;
       const auto passes = 2 + static_cast<int>(auto_gamma); // log-Gaussian, -f', and the gamma scan
       report << "Performance: wall=" << wall_seconds << "s";
       if (cpu_seconds && wall_seconds > 0.0)
         report << " cpu=" << *cpu_seconds << "s effective_parallelism=" << *cpu_seconds / wall_seconds << "x";
       else
         report << " cpu=n/a effective_parallelism=n/a";
       report << " throughput=" << throughput << " points/s workers=" << std::min(jobs, mesh.size())
              << " passes=" << passes << " cells=" << nr_cells << " groups=" << positive.size() + negative.size() << '\n';
       report << "Phases:";
       for (const auto &[name, seconds] : phases) report << " " << name << "=" << seconds << "s";
       report << '\n';
     }
     std::cout << report.str();
     NRG::Tools::finish_output(std::cout, "<stdout>");
   }

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
       << " -X -- choose gamma automatically (Sec. V.2): scan gamma/T in [0.01, 1], write gamma_scan.dat, take the\n"
       << "      smallest gamma at which eps_rms falls below the threshold\n"
       << " -e <threshold> -- threshold on eps_rms for -X (default 0.002)\n"
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
     while ((c_ = getopt_long(argc, argv, "vj:m:M:r:L:oIx:S:Xe:", long_options, nullptr)) != -1) {
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
       case 'X': auto_gamma = true; break;
       case 'e': eps_threshold = std::atof(optarg); break;
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
     if (!(eps_threshold > 0.0) || !std::isfinite(eps_threshold)) throw std::invalid_argument("The -e threshold must be finite and greater than 0.");
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
     std::map<int, std::vector<Cell>> by_width_pos, by_width_neg; // key: floor(log2(width))
     for (const auto &[key, w] : cells) {
       const auto &[omega, sigma] = key;
       const auto weight = w / Nz;
       total_weight += weight;
       if (omega == 0.0 || weight == 0.0) continue;
       const auto width = alpha * std::max(sigma, sigma_min) / lnL;
       if (!(width > 0.0) || !std::isfinite(width)) throw std::runtime_error("A cell has a non-positive broadening width; use -S.");
       const auto k = static_cast<int>(std::floor(std::log2(width)));
       (omega > 0.0 ? by_width_pos : by_width_neg)[k].push_back({std::log(std::abs(omega)), omega, width, weight});
     }
     const auto by_log = [](const Cell &x, const Cell &y) { return x.logabs < y.logabs; };
     const auto make_groups = [&by_log](std::map<int, std::vector<Cell>> &by_width, std::vector<Group> &groups) {
       for (auto &[k, group_cells] : by_width) {
         std::sort(group_cells.begin(), group_cells.end(), by_log);
         double b = 0.0;
         for (const auto &c : group_cells) b = std::max(b, c.width);
         groups.push_back({b * (7.0 + b / 4.0), std::move(group_cells)});
       }
     };
     make_groups(by_width_pos, positive);
     make_groups(by_width_neg, negative);
     nr_cells = cells.size();
     std::cout << cells.size() << " cells, total weight=" << total_weight << std::endl;
     if (verbose) std::cerr << positive.size() + negative.size() << " groups of cells by kernel width" << std::endl;
   }

   // Each cell broadened with its own log-Gaussian (Eq. 17a). Cells further than b (7 + b/4) in ln|omega| from the
   // output point, b the widest kernel of their group, contribute less than exp(-49) relative and are skipped.
   std::vector<double> broaden_log() const {
     std::vector<double> result(mesh.size());
     for_each_index(mesh.size(), jobs, [&](const size_t index) {
       const auto m = mesh[index];
       const auto &groups = m > 0.0 ? positive : negative;
       const auto logm = std::log(std::abs(m));
       double sum = 0.0;
       for (const auto &group : groups) {
         const auto &cells = group.cells;
         const auto lo = std::lower_bound(cells.begin(), cells.end(), logm - group.window,
                                          [](const Cell &c, const double v) { return c.logabs < v; });
         for (auto it = lo; it != cells.end() && it->logabs <= logm + group.window; ++it)
           sum += it->weight * NRG::Broadening::log_gaussian(m, it->omega, it->width);
       }
       result[index] = sum;
     });
     return result;
   }

   // The -f' convolution of b at a single frequency x, with the quadrature of convolve() (broaden_common.hpp).
   [[nodiscard]] double convolve_at(const std::vector<double> &b, const double x, const double sigma) const {
     double sum = 0.0;
     for (size_t j = 1; j + 1 < mesh.size(); j++) {
       const auto width = (mesh[j + 1] - mesh[j - 1]) / 2.0;
       sum += b[j] * derfd_kernel(x, mesh[j], sigma) * width;
     }
     return sum;
   }

   // eps_rms of Eq. 24: the normalized standard error of a least-squares quadratic fit through the points (x_i, y_i).
   template<size_t n>
   [[nodiscard]] static double eps_rms(const std::array<double, n> &x, const std::array<double, n> &y) {
     static_assert(n > 3);
     const auto scale = std::max(std::abs(x.front()), std::abs(x.back())); // fit in u = x/scale, for conditioning
     std::array<double, 5> su{};                                          // sums of u^k, k = 0..4
     std::array<double, 3> suy{};                                         // sums of u^k y, k = 0..2
     for (size_t i = 0; i < n; i++) {
       const auto u = x[i] / scale;
       double uk = 1.0;
       for (size_t k = 0; k < 5; k++) {
         su[k] += uk;
         if (k < 3) suy[k] += uk * y[i];
         uk *= u;
       }
     }
     // Normal equations sum_l su[k+l] c_l = suy[k], solved by Gaussian elimination with partial pivoting.
     std::array<std::array<double, 4>, 3> m{};
     for (size_t k = 0; k < 3; k++) {
       for (size_t l = 0; l < 3; l++) m[k][l] = su[k + l];
       m[k][3] = suy[k];
     }
     for (size_t col = 0; col < 3; col++) {
       size_t pivot = col;
       for (size_t r = col + 1; r < 3; r++)
         if (std::abs(m[r][col]) > std::abs(m[pivot][col])) pivot = r;
       std::swap(m[col], m[pivot]);
       for (size_t r = col + 1; r < 3; r++) {
         const auto f = m[r][col] / m[col][col];
         for (size_t c = col; c < 4; c++) m[r][c] -= f * m[col][c];
       }
     }
     std::array<double, 3> c{};
     for (size_t k = 3; k-- > 0;) {
       double v = m[k][3];
       for (size_t l = k + 1; l < 3; l++) v -= m[k][l] * c[l];
       c[k] = v / m[k][k];
     }
     double residual = 0.0, mean = 0.0;
     for (size_t i = 0; i < n; i++) {
       const auto u = x[i] / scale;
       const auto fit = c[0] + c[1] * u + c[2] * u * u;
       residual += (y[i] - fit) * (y[i] - fit);
       mean += y[i];
     }
     mean /= n;
     if (!(mean > 0.0))
       throw std::runtime_error("The spectrum near omega=0 is not positive, so eps_rms is undefined; set gamma with -x.");
     return std::sqrt(residual / (n - 3)) / mean;
   }

   // Automatic gamma (Sec. V.2 of the paper): for gamma/T on a logarithmic grid over [0.01, 1], the -f' convolution of
   // the log-broadened curve b at m_l = 11 equally spaced points on [-gamma Lambda^2, gamma Lambda^2] (n_l = 2), a
   // quadratic fit and eps_rms. The result is the smallest gamma at which eps_rms falls below the threshold, coming
   // from larger values, interpolated linearly in (ln gamma, ln eps_rms) between the two grid points around the
   // crossing; if the first point is already below, that point; if no point is, the minimum of eps_rms, with a
   // warning. Rule (iv) of the paper (the geometric mean across an irregular window) is left to the user, who sees
   // the scan in gamma_scan.dat and can set gamma with -x.
   [[nodiscard]] double scan_gamma(const std::vector<double> &b) const {
     constexpr size_t npoints = 41;
     constexpr size_t ml = 11;
     const auto mesh_min = std::abs(*std::min_element(mesh.begin(), mesh.end(), [](const double x, const double y) { return std::abs(x) < std::abs(y); }));
     if (mesh_min > 0.01 * T)
       throw std::runtime_error("-X: the output mesh must reach below 0.01*T (use -m).");
     std::vector<double> ratios(npoints), eps(npoints);
     for (size_t k = 0; k < npoints; k++) {
       ratios[k] = 0.01 * std::pow(100.0, double(k) / double(npoints - 1));
       const auto gamma = ratios[k] * T;
       const auto omega_l = gamma * Lambda * Lambda;
       std::array<double, ml> x{}, y{};
       for (size_t i = 0; i < ml; i++) {
         x[i] = omega_l * (double(i) / 5.0 - 1.0);
         y[i] = convolve_at(b, x[i], gamma);
       }
       eps[k] = eps_rms(x, y);
     }
     {
       std::ofstream F(scan_filename);
       if (!F) throw std::runtime_error("Failed to open " + scan_filename + " for writing.");
       F << std::setprecision(17) << "# gamma/T eps_rms\n";
       for (size_t k = 0; k < npoints; k++) F << ratios[k] << " " << eps[k] << '\n';
       NRG::Tools::finish_output(F, scan_filename);
     }
     if (eps[0] <= eps_threshold) return ratios[0];
     for (size_t k = 1; k < npoints; k++) {
       if (eps[k] <= eps_threshold) {
         if (!(eps[k] > 0.0)) return ratios[k];
         const auto l0 = std::log(ratios[k - 1]), l1 = std::log(ratios[k]);
         const auto e0 = std::log(eps[k - 1]), e1 = std::log(eps[k]);
         return std::exp(l0 + (std::log(eps_threshold) - e0) * (l1 - l0) / (e1 - e0));
       }
     }
     const auto best = std::min_element(eps.begin(), eps.end()) - eps.begin();
     std::cerr << "abroaden: warning: eps_rms never falls below " << eps_threshold << "; taking the minimum, "
               << eps[best] << " at gamma/T=" << ratios[best] << ". See " << scan_filename << "." << std::endl;
     return ratios[best];
   }

 public:
   AdaptiveBroaden(int argc, char *argv[]) : wall_start{std::chrono::steady_clock::now()}, cpu_start{std::clock()} {
     std::cout << "abroaden - adaptive broadening tool" << std::endl;
     cmd_line(argc, argv);
   }
   void calc() {
     timed("read", [this] { read_and_merge(); });
     mesh = mesh_filename.empty() ? make_mesh(broaden_min, broaden_max, broaden_ratio, 0.0, true, true)
                                  : load_mesh(mesh_filename, verbose);
     if (mesh.empty()) throw std::runtime_error("Output mesh is empty.");
     if (verbose)
       std::cerr << "alpha=" << alpha << " T=" << T << " Lambda=" << Lambda << " gamma=" << gamma_ratio * T
                 << " sigma_min=" << sigma_min << " mesh points=" << mesh.size() << " jobs=" << jobs << std::endl;
     a = timed("log_gaussian", [this] { return broaden_log(); });
     if (auto_gamma) {
       gamma_ratio = timed("gamma_scan", [this] { return scan_gamma(a); });
       std::cout << std::setprecision(17) << "Automatic gamma/T=" << gamma_ratio << " (eps_rms threshold " << eps_threshold
                 << ", scan in " << scan_filename << ")" << std::endl;
     }
     timed("derfd", [this] { convolve(mesh, a, gamma_ratio * T, jobs, derfd_kernel); }); // Eq. 16: after the log-Gaussian, over all frequencies
     std::cout << "Estimated weight (trapezoidal rule)=" << trapez(mesh, a) << std::endl;
     save(output_filename, mesh, a, verbose);
     report_timing();
   }
};

} // namespace NRG::Broaden

#endif
