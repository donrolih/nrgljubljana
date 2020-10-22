#ifndef _nrg_lib_h_
#define _nrg_lib_h_

// This header is included for both executable and library.

#define NRG_COMMON

#include <string>
#include "workdir.h"

void run_nrg_master(const Workdir &workdir, const bool embedded);
void run_nrg_slave(); // note: only defined if compled using NRG_MPI
void print_about_message();

#endif
