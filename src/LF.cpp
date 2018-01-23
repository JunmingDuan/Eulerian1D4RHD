#include "Eulerian1D.h"

bU Eulerian1D::LF(const bU& CONL, const bU& CONR, const bU& PRIL, const bU& PRIR, const double alpha) {
  return 0.5*(F(CONL, PRIL) + F(CONR, PRIR)) - 0.5*(CONR-CONL)*alpha;
}

void Eulerian1D::cal_flux_LF(Sol& ReconL_Con, Sol& ReconR_Con, Sol& ReconL_Pri, Sol& ReconR_Pri,
    Sol& FLUX, double alpha) {
#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x+1; ++i) {
    FLUX[i] = LF(ReconL_Con[i], ReconR_Con[i], ReconL_Pri[i], ReconR_Pri[i], alpha);
  }
}

