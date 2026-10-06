#include <cstdlib>
#include <exception>
#include <iostream>
#include <ostream>

#include <common/version.hpp>

#include "abroaden.hpp"

int main(int argc, char *argv[]) {
  if (NRG::Tools::report_version_if_requested(argc, argv, "abroaden")) return EXIT_SUCCESS;
  try {
    NRG::Broaden::AdaptiveBroaden abroaden(argc, argv);
    abroaden.calc();
  } catch (const std::exception &e) {
    std::cerr << "abroaden: error: " << e.what() << std::endl;
    return 1;
  } catch (...) {
    std::cerr << "abroaden: error: unknown exception" << std::endl;
    return 1;
  }
}
