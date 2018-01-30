#include "Eulerian1D.h"
#include "WENO_nonuniform.h"
#include "ENO.h"

#define SECOND_ORDER_ENO
//#define THIRD_ORDER_WENO
#ifdef SECOND_ORDER_ENO
void Eulerian1D::Reconstruction(const Sol& sol, const VEC& mesh,
    Sol& ReconL_Con, Sol& ReconR_Con, Sol& ReconL_Pri, Sol& ReconR_Pri) { //ENO2
  if(is_RECON == 0) {//sol,gl,gr are conservative variables
#pragma omp parallel for num_threads(Nthread)
    for(u_int i = 0; i < N_x+1; ++i) {
      if(i == 0) {
        if(BD == 1) {
          ReconL_Con[i] = sol[0];
        }
        else if(BD == 2) {
          ReconL_Con[i] = sol[N_x-1];
        }
        ReconR_Con[i] = sol[i];
        ReconL_Pri[i] = Con2Pri(ReconL_Con[i], Gamma[N_x-1]);
        ReconR_Pri[i] = Con2Pri(ReconR_Con[i], Gamma[i]);
      }
      else if(i == N_x) {
        if(BD == 1) {
          ReconR_Con[i] = sol[N_x-1];
        }
        else if(BD == 2) {
          ReconR_Con[i] = sol[0];
        }
        ReconL_Con[i] = sol[i-1];
        ReconL_Pri[i] = Con2Pri(ReconL_Con[i], Gamma[i-1]);
        ReconR_Pri[i] = Con2Pri(ReconR_Con[i], Gamma[0]);
      }
      else {
        ReconL_Con[i] = sol[i-1];
        ReconR_Con[i] = sol[i];
        ReconL_Pri[i] = Con2Pri(ReconL_Con[i], Gamma[i-1]);
        ReconR_Pri[i] = Con2Pri(ReconR_Con[i], Gamma[i]);
      }
    }

  }
  else if(is_RECON == 1) {//sol,gl,gr are conservative variables
    VEC h(N_x);
#pragma omp parallel for num_threads(Nthread)
    for(u_int i = 0; i < N_x; ++i) {
      h[i] = mesh[i+1] - mesh[i];
    }
#pragma omp parallel for num_threads(Nthread)
    for(u_int i = 0; i < N_x; ++i) {//number of the cell
      for(u_int d = 0; d < 3; ++d) {
        if(i == 0) {
          if(BD == 1) {
            ENO2(h[0], h[i], h[i+1], sol[0][d], sol[i][d], sol[i+1][d], ReconR_Con[i][d], ReconL_Con[i+1][d]);
          }
          else if(BD == 2) {
            ENO2(h[N_x-1], h[i], h[i+1], sol[N_x-1][d], sol[i][d], sol[i+1][d], ReconR_Con[i][d], ReconL_Con[i+1][d]);
          }
        }
        else if(i == N_x-1) {
          if(BD == 1) {
            ENO2(h[i-1], h[i], h[N_x-1], sol[i-1][d], sol[i][d], sol[N_x-1][d], ReconR_Con[i][d], ReconL_Con[i+1][d]);
          }
          else if(BD == 2) {
            ENO2(h[i-1], h[i], h[0], sol[i-1][d], sol[i][d], sol[0][d], ReconR_Con[i][d], ReconL_Con[i+1][d]);
          }
        }
        else {
          ENO2(h[i-1], h[i], h[i+1], sol[i-1][d], sol[i][d], sol[i+1][d], ReconR_Con[i][d], ReconL_Con[i+1][d]);
        }
      }
      ReconL_Pri[i+1] = Con2Pri(ReconL_Con[i+1], Gamma[i]);
      ReconR_Pri[i] = Con2Pri(ReconR_Con[i], Gamma[i]);
    }
    if(BD == 1) {//outflow
      ReconL_Con[0] = ReconR_Con[0];
      ReconR_Con[N_x] = ReconL_Con[N_x];
      ReconL_Pri[0] = ReconR_Pri[0];
      ReconR_Pri[N_x] = ReconL_Pri[N_x];
    }
    else if(BD == 2) {//period
      ReconL_Con[0] = ReconL_Con[N_x];
      ReconR_Con[N_x] = ReconR_Con[0];
      ReconL_Pri[0] = ReconL_Pri[N_x];
      ReconR_Pri[N_x] = ReconR_Pri[0];
    }
  }
  else if(is_RECON == 2) {//sol is primitive variables
    VEC h(N_x);
#pragma omp parallel for num_threads(Nthread)
    for(u_int i = 0; i < N_x; ++i) {
      h[i] = mesh[i+1] - mesh[i];
    }
#pragma omp parallel for num_threads(Nthread)
    for(u_int i = 0; i < N_x; ++i) {//number of the cell
      for(u_int d = 0; d < 3; ++d) {
        if(i == 0) {
          if(BD == 1) {
            ENO2(h[0], h[i], h[i+1], sol[0][d], sol[i][d], sol[i+1][d], ReconR_Pri[i][d], ReconL_Pri[i+1][d]);
          }
          else if(BD == 2) {
            ENO2(h[N_x-1], h[i], h[i+1], sol[N_x-1][d], sol[i][d], sol[i+1][d], ReconR_Pri[i][d], ReconL_Pri[i+1][d]);
          }
        }
        else if(i == N_x-1) {
          if(BD == 1) {
            ENO2(h[i-1], h[i], h[N_x-1], sol[i-1][d], sol[i][d], sol[N_x-1][d], ReconR_Pri[i][d], ReconL_Pri[i+1][d]);
          }
          else if(BD == 2) {
            ENO2(h[i-1], h[i], h[0], sol[i-1][d], sol[i][d], sol[0][d], ReconR_Pri[i][d], ReconL_Pri[i+1][d]);
          }
        }
        else {
          ENO2(h[i-1], h[i], h[i+1], sol[i-1][d], sol[i][d], sol[i+1][d], ReconR_Pri[i][d], ReconL_Pri[i+1][d]);
        }
      }
      ReconL_Con[i+1] = Pri2Con(ReconL_Pri[i+1], Gamma[i]);
      ReconR_Con[i] = Pri2Con(ReconR_Pri[i], Gamma[i]);
    }
    if(BD == 1) {//outflow
      ReconL_Con[0] = ReconR_Con[0];
      ReconR_Con[N_x] = ReconL_Con[N_x];
      ReconL_Pri[0] = ReconR_Pri[0];
      ReconR_Pri[N_x] = ReconL_Pri[N_x];
    }
    else if(BD == 2) {//period
      ReconL_Con[0] = ReconL_Con[N_x];
      ReconR_Con[N_x] = ReconR_Con[0];
      ReconL_Pri[0] = ReconL_Pri[N_x];
      ReconR_Pri[N_x] = ReconR_Pri[0];
    }
  }
  else if(is_RECON == 3) {//sol is primitive variables
    VEC h(N_x);
    Sol CON_tmp(N_x);
    bU U1, U2, U3, U4;
    double h1, h2, h3, h4;
    MAT LMAT, RMAT;
#pragma omp parallel for num_threads(Nthread)
    for(u_int i = 0; i < N_x; ++i) {
      h[i] = mesh[i+1] - mesh[i];
      CON_tmp[i] = Pri2Con(sol[i], Gamma[i]);
    }
#pragma omp parallel for num_threads(Nthread)
    for(u_int i = 1; i < N_x; ++i) {//perform characteristic decomposition
      if(i == 1) {
        ROE_AV_MAT(sol[i-1], sol[i], Gamma[i-1], Gamma[i], RMAT, LMAT);
        h1 = h[0];
        h2 = h[i-1];
        h3 = h[i];
        h4 = h[i+1];
        if(BD == 1) {
          U1 = multiply(CON_tmp[0], LMAT);
        }
        else if(BD == 2) {
          U1 = multiply(CON_tmp[N_x-1], LMAT);
        }
        U2 = multiply(CON_tmp[i-1], LMAT);
        U3 = multiply(CON_tmp[i], LMAT);
        U4 = multiply(CON_tmp[i+1], LMAT);
        for(u_int d = 0; d < 3; ++d) {
          ReconR_Con[i][d] = ENO2_CELL_L(h2, h3, h4, U2[d], U3[d], U4[d]);
          ReconL_Con[i][d] = ENO2_CELL_R(h1, h2, h3, U1[d], U2[d], U3[d]);
        }
        ReconL_Con[i] = multiply(ReconL_Con[i], RMAT);
        ReconR_Con[i] = multiply(ReconR_Con[i], RMAT);
        ReconL_Pri[i] = Con2Pri(ReconL_Con[i], Gamma[i-1]);
        ReconR_Pri[i] = Con2Pri(ReconR_Con[i], Gamma[i]);
      }
      else if(i == N_x-1) {
        ROE_AV_MAT(sol[i-1], sol[i], Gamma[i-1], Gamma[i], RMAT, LMAT);
        h1 = h[i-2];
        h2 = h[i-1];
        h3 = h[i];
        h4 = h[N_x-1];
        U1 = multiply(CON_tmp[i-2], LMAT);
        U2 = multiply(CON_tmp[i-1], LMAT);
        U3 = multiply(CON_tmp[i], LMAT);
        if(BD == 1) {
          U4 = multiply(CON_tmp[N_x-1], LMAT);
        }
        else if(BD == 2) {
          U4 = multiply(CON_tmp[0], LMAT);
        }
        for(u_int d = 0; d < 3; ++d) {
          ReconR_Con[i][d] = ENO2_CELL_L(h2, h3, h4, U2[d], U3[d], U4[d]);
          ReconL_Con[i][d] = ENO2_CELL_R(h1, h2, h3, U1[d], U2[d], U3[d]);
        }
        ReconL_Con[i] = multiply(ReconL_Con[i], RMAT);
        ReconR_Con[i] = multiply(ReconR_Con[i], RMAT);
        ReconL_Pri[i] = Con2Pri(ReconL_Con[i], Gamma[i-1]);
        ReconR_Pri[i] = Con2Pri(ReconR_Con[i], Gamma[i]);
      }
      else {
        ROE_AV_MAT(sol[i-1], sol[i], Gamma[i-1], Gamma[i], RMAT, LMAT);
        h1 = h[i-2];
        h2 = h[i-1];
        h3 = h[i];
        h4 = h[i+1];
        U1 = multiply(CON_tmp[i-2], LMAT);
        U2 = multiply(CON_tmp[i-1], LMAT);
        U3 = multiply(CON_tmp[i], LMAT);
        U4 = multiply(CON_tmp[i+1], LMAT);
        for(u_int d = 0; d < 3; ++d) {
          ReconR_Con[i][d] = ENO2_CELL_L(h2, h3, h4, U2[d], U3[d], U4[d]);
          ReconL_Con[i][d] = ENO2_CELL_R(h1, h2, h3, U1[d], U2[d], U3[d]);
        }
        ReconL_Con[i] = multiply(ReconL_Con[i], RMAT);
        ReconR_Con[i] = multiply(ReconR_Con[i], RMAT);
        ReconL_Pri[i] = Con2Pri(ReconL_Con[i], Gamma[i-1]);
        ReconR_Pri[i] = Con2Pri(ReconR_Con[i], Gamma[i]);
      }
   }
    //deal with left BD
    {
      if(BD == 1) {
        ROE_AV_MAT(sol[0], sol[0], Gamma[0], Gamma[0], RMAT, LMAT);
      }
      else if(BD == 2) {
        ROE_AV_MAT(sol[N_x-1], sol[0], Gamma[N_x-1], Gamma[0], RMAT, LMAT);
      }
      h2 = h[0];
      h3 = h[0];
      h4 = h[1];
      U2 = multiply(CON_tmp[0], LMAT);
      U3 = multiply(CON_tmp[0], LMAT);
      U4 = multiply(CON_tmp[1], LMAT);
      for(u_int d = 0; d < 3; ++d) {
        ReconR_Con[0][d] = ENO2_CELL_L(h2, h3, h4, U2[d], U3[d], U4[d]);
      }
      ReconR_Con[0] = multiply(ReconR_Con[0], RMAT);
      ReconR_Pri[0] = Con2Pri(ReconR_Con[0], Gamma[0]);
    }
    //deal with right BD
    {
      if(BD == 1) {
        ROE_AV_MAT(sol[N_x-1], sol[N_x-1], Gamma[N_x-1], Gamma[N_x-1], RMAT, LMAT);
      }
      else if(BD == 2) {
        ROE_AV_MAT(sol[N_x-1], sol[0], Gamma[N_x-1], Gamma[0], RMAT, LMAT);
      }
      h1 = h[N_x-2];
      h2 = h[N_x-1];
      h3 = h[N_x-1];
      U1 = multiply(CON_tmp[N_x-2], LMAT);
      U2 = multiply(CON_tmp[N_x-1], LMAT);
      U3 = multiply(CON_tmp[N_x-1], LMAT);
      for(u_int d = 0; d < 3; ++d) {
        ReconL_Con[N_x][d] = ENO2_CELL_R(h1, h2, h3, U1[d], U2[d], U3[d]);
      }
      ReconL_Con[N_x] = multiply(ReconL_Con[N_x], RMAT);
      ReconL_Pri[N_x] = Con2Pri(ReconL_Con[N_x], Gamma[N_x-1]);
    }
    if(BD == 1) {//outflow
      ReconL_Con[0] = ReconR_Con[0];
      ReconR_Con[N_x] = ReconL_Con[N_x];
      ReconL_Pri[0] = ReconR_Pri[0];
      ReconR_Pri[N_x] = ReconL_Pri[N_x];
    }
    else if(BD == 2) {//period
      ReconL_Con[0] = ReconL_Con[N_x];
      ReconR_Con[N_x] = ReconR_Con[0];
      ReconL_Pri[0] = ReconL_Pri[N_x];
      ReconR_Pri[N_x] = ReconR_Pri[0];
    }
  }

}
#endif

#ifdef THIRD_ORDER_WENO
void Eulerian1D::Reconstruction(const Sol& sol, const VEC& mesh,
    Sol& ReconL_Con, Sol& ReconR_Con, Sol& ReconL_Pri, Sol& ReconR_Pri) { //WENO
  if(is_RECON == 0) {//sol,gl,gr are conservative variables
#pragma omp parallel for num_threads(Nthread)
    for(u_int i = 0; i < N_x+1; ++i) {
      if(i == 0) {
        if(BD == 1) {
          ReconL_Con[i] = sol[0];
        }
        else if(BD == 2) {
          ReconL_Con[i] = sol[N_x-1];
        }
        ReconR_Con[i] = sol[i];
        ReconL_Pri[i] = Con2Pri(ReconL_Con[i], Gamma[N_x-1]);
        ReconR_Pri[i] = Con2Pri(ReconR_Con[i], Gamma[i]);
      }
      else if(i == N_x) {
        if(BD == 1) {
          ReconR_Con[i] = sol[N_x-1];
        }
        else if(BD == 2) {
          ReconR_Con[i] = sol[0];
        }
        ReconL_Con[i] = sol[i-1];
        ReconL_Pri[i] = Con2Pri(ReconL_Con[i], Gamma[i-1]);
        ReconR_Pri[i] = Con2Pri(ReconR_Con[i], Gamma[0]);
      }
      else {
        ReconL_Con[i] = sol[i-1];
        ReconR_Con[i] = sol[i];
        ReconL_Pri[i] = Con2Pri(ReconL_Con[i], Gamma[i-1]);
        ReconR_Pri[i] = Con2Pri(ReconR_Con[i], Gamma[i]);
      }
    }
  }
  else if(is_RECON == 1) {//sol,gl,gr are conservative variables
    VEC h(N_x);
#pragma omp parallel for num_threads(Nthread)
    for(u_int i = 0; i < N_x; ++i) {
      h[i] = mesh[i+1] - mesh[i];
    }
#pragma omp parallel for num_threads(Nthread)
    for(u_int i = 0; i < N_x; ++i) {//number of the cell
      for(u_int d = 0; d < 3; ++d) {
        if(i == 0) {
          if(BD == 1) {
            WENO3(h[0], h[i], h[i+1], sol[0][d], sol[i][d], sol[i+1][d], ReconR_Con[i][d], ReconL_Con[i+1][d]);
          }
          else if(BD == 2) {
            WENO3(h[N_x-1], h[i], h[i+1], sol[N_x-1][d], sol[i][d], sol[i+1][d], ReconR_Con[i][d], ReconL_Con[i+1][d]);
          }
        }
        else if(i == N_x-1) {
          if(BD == 1) {
            WENO3(h[i-1], h[i], h[N_x-1], sol[i-1][d], sol[i][d], sol[N_x-1][d], ReconR_Con[i][d], ReconL_Con[i+1][d]);
          }
          else if(BD == 2) {
            WENO3(h[i-1], h[i], h[0], sol[i-1][d], sol[i][d], sol[0][d], ReconR_Con[i][d], ReconL_Con[i+1][d]);
          }
        }
        else {
          WENO3(h[i-1], h[i], h[i+1], sol[i-1][d], sol[i][d], sol[i+1][d], ReconR_Con[i][d], ReconL_Con[i+1][d]);
        }
      }
      ReconL_Pri[i+1] = Con2Pri(ReconL_Con[i+1], Gamma[i]);
      ReconR_Pri[i] = Con2Pri(ReconR_Con[i], Gamma[i]);
    }
    if(BD == 1) {//outflow
      ReconL_Con[0] = ReconL_Con[0];
      ReconR_Con[N_x] = ReconR_Con[N_x];
      ReconL_Pri[0] = ReconL_Pri[0];
      ReconR_Pri[N_x] = ReconR_Pri[N_x];
    }
    else if(BD == 2) {//period
      ReconL_Con[0] = ReconL_Con[N_x];
      ReconR_Con[N_x] = ReconR_Con[0];
      ReconL_Pri[0] = ReconL_Pri[N_x];
      ReconR_Pri[N_x] = ReconR_Pri[0];
    }
  }
  else if(is_RECON == 2) {//sol is primitive variables
    VEC h(N_x);
#pragma omp parallel for num_threads(Nthread)
    for(u_int i = 0; i < N_x; ++i) {
      h[i] = mesh[i+1] - mesh[i];
    }
#pragma omp parallel for num_threads(Nthread)
    for(u_int i = 0; i < N_x; ++i) {//number of the cell
      for(u_int d = 0; d < 3; ++d) {
        if(i == 0) {
          if(BD == 1) {
            WENO3(h[0], h[i], h[i+1], sol[0][d], sol[i][d], sol[i+1][d], ReconR_Pri[i][d], ReconL_Pri[i+1][d]);
          }
          else if(BD == 2) {
            WENO3(h[N_x-1], h[i], h[i+1], sol[N_x-1][d], sol[i][d], sol[i+1][d], ReconR_Pri[i][d], ReconL_Pri[i+1][d]);
          }
        }
        else if(i == N_x-1) {
          if(BD == 1) {
            WENO3(h[i-1], h[i], h[N_x-1], sol[i-1][d], sol[i][d], sol[N_x-1][d], ReconR_Pri[i][d], ReconL_Pri[i+1][d]);
          }
          else if(BD == 2) {
            WENO3(h[i-1], h[i], h[0], sol[i-1][d], sol[i][d], sol[0][d], ReconR_Pri[i][d], ReconL_Pri[i+1][d]);
          }
        }
        else {
          WENO3(h[i-1], h[i], h[i+1], sol[i-1][d], sol[i][d], sol[i+1][d], ReconR_Pri[i][d], ReconL_Pri[i+1][d]);
        }
      }
      ReconL_Con[i+1] = Pri2Con(ReconL_Pri[i+1], Gamma[i]);
      ReconR_Con[i] = Pri2Con(ReconR_Pri[i], Gamma[i]);
    }
    if(BD == 1) {//outflow
      ReconL_Con[0] = ReconR_Con[0];
      ReconR_Con[N_x] = ReconL_Con[N_x];
      ReconL_Pri[0] = ReconR_Pri[0];
      ReconR_Pri[N_x] = ReconL_Pri[N_x];
    }
    else if(BD == 2) {//period
      ReconL_Con[0] = ReconL_Con[N_x];
      ReconR_Con[N_x] = ReconR_Con[0];
      ReconL_Pri[0] = ReconL_Pri[N_x];
      ReconR_Pri[N_x] = ReconR_Pri[0];
    }
  }

}
#endif

