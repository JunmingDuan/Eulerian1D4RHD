#include "Eulerian1D.h"

bU Eulerian1D::LLF(const bU& CONL, const bU& CONR, const bU& PRIL, const bU& PRIR, const double alpha) {
  return 0.5*(F(CONL, PRIL) + F(CONR, PRIR)) - 0.5*(CONR-CONL)*alpha;
}

void Eulerian1D::cal_flux_LLF(Sol& ReconL_Con, Sol& ReconR_Con, Sol& ReconL_Pri, Sol& ReconR_Pri, Sol& FLUX) {
#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x+1; ++i) {
    double alpha, laml, lamr;
    //cal local characteristic speed
    if(i == 0) {
      laml = cal_max_lambda_Eul(ReconL_Con[i], ReconL_Pri[i], Gamma[0]);
      lamr = cal_max_lambda_Eul(ReconR_Con[i], ReconR_Pri[i], Gamma[i]);
    }
    else if(i == N_x) {
      laml = cal_max_lambda_Eul(ReconL_Con[i], ReconL_Pri[i], Gamma[i-1]);
      lamr = cal_max_lambda_Eul(ReconR_Con[i], ReconR_Pri[i], Gamma[N_x-1]);
    }
    else {
      laml = cal_max_lambda_Eul(ReconL_Con[i], ReconL_Pri[i], Gamma[i-1]);
      lamr = cal_max_lambda_Eul(ReconR_Con[i], ReconR_Pri[i], Gamma[i]);
    }
    alpha = std::max(laml, lamr);
    //alpha *= 1.3;
    FLUX[i] = LLF(ReconL_Con[i], ReconR_Con[i], ReconL_Pri[i], ReconR_Pri[i], alpha);
  }
}

