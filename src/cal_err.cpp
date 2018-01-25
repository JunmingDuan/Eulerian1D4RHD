#include "Eulerian1D.h"

double Eulerian1D::cal_err(int l) {
  Sol exact(N_x);
  double err(0);
  for(u_int i = 0; i < N_x; ++i) {
    exact[i] = initial(t_end, 0.5*(mesh[i]+mesh[i+1]), Gamma[i]);
  }
  if(l == 1) {
    for(u_int i = 0; i < N_x; ++i) {
      err += fabs(Pri[i][0]-exact[i][0]);
    }
    return err/N_x;
  }
  else if(l == 2) {
    for(u_int i = 0; i < N_x; ++i) {
      err += pow(Pri[i][0]-exact[i][0], 2);
    }
    return sqrt(err)/N_x;
  }
  else if(l == -1) {
    for(u_int i = 0; i < N_x; ++i) {
      err = std::max(err, fabs(Pri[i][0]-exact[i][0]));
    }
    return err;
  }
  return 0;
}

