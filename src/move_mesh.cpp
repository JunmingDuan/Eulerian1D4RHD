#include "Eulerian1D.h"

void Eulerian1D::cal_us_roeav(Sol& ReconL_Pri, Sol& ReconR_Pri, vvector<double>& us) {
#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x+1; ++i) {
    us[i] = 0;
  }
}

void Eulerian1D::move_mesh(vvector<double>& mesh, vvector<double>& us, double dt, vvector<double>& mesh1) {
#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x+1; ++i) {
    mesh1[i] = mesh[i] + dt * us[i];
  }
}

