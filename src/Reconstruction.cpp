#include "Eulerian1D.h"
#include "WENO_nonuniform.h"
#include "ENO.h"

//ENO
//void Eulerian1D::Reconstruction(const Sol& sol, const VEC& mesh,
    //Sol& ReconL_Con, Sol& ReconR_Con, Sol& ReconL_Pri, Sol& ReconR_Pri) {
  //if(is_RECON == 0) {//sol,gl,gr are conservative variables
//#pragma omp parallel for num_threads(Nthread)
    //for(u_int i = 0; i < N_x+1; ++i) {
      //if(i == 0) {
        //if(BD == 2) {
          //ReconL_Con[i] = sol[N_x-1];
        //}
        //ReconR_Con[i] = sol[i];
        //ReconL_Pri[i] = Con2Pri(ReconL_Con[i], Gamma[N_x-1]);
        //ReconR_Pri[i] = Con2Pri(ReconR_Con[i], Gamma[i]);
      //}
      //else if(i == N_x) {
        //if(BD == 2) {
          //ReconR_Con[i] = sol[0];
        //}
        //ReconL_Con[i] = sol[i-1];
        //ReconL_Pri[i] = Con2Pri(ReconL_Con[i], Gamma[i-1]);
        //ReconR_Pri[i] = Con2Pri(ReconR_Con[i], Gamma[0]);
      //}
      //else {
        //ReconL_Con[i] = sol[i-1];
        //ReconR_Con[i] = sol[i];
        //ReconL_Pri[i] = Con2Pri(ReconL_Con[i], Gamma[i-1]);
        //ReconR_Pri[i] = Con2Pri(ReconR_Con[i], Gamma[i]);
      //}
    //}
  //}
  //else if(is_RECON == 1) {//sol,gl,gr are conservative variables
    //VEC h(N_x);
//#pragma omp parallel for num_threads(Nthread)
    //for(u_int i = 0; i < N_x; ++i) {
      //h[i] = mesh[i+1] - mesh[i];
    //}
//#pragma omp parallel for num_threads(Nthread)
    //for(u_int i = 0; i < N_x; ++i) {//number of the cell
      //for(u_int d = 0; d < 3; ++d) {
        //if(i == 0) {
          //if(BD == 2) {
            //ENO2(h[N_x-1], h[i], h[i+1], sol[N_x-1][d], sol[i][d], sol[i+1][d], ReconR_Con[i][d], ReconL_Con[i+1][d]);
          //}
        //}
        //else if(i == N_x-1) {
          //if(BD == 2) {
            //ENO2(h[i-1], h[i], h[0], sol[i-1][d], sol[i][d], sol[0][d], ReconR_Con[i][d], ReconL_Con[i+1][d]);
          //}
        //}
        //else {
          //ENO2(h[i-1], h[i], h[i+1], sol[i-1][d], sol[i][d], sol[i+1][d], ReconR_Con[i][d], ReconL_Con[i+1][d]);
        //}
      //}
      //ReconL_Pri[i+1] = Con2Pri(ReconL_Con[i+1], Gamma[i]);
      //ReconR_Pri[i] = Con2Pri(ReconR_Con[i], Gamma[i]);
    //}
    //if(BD == 2) {//period
      //ReconL_Con[0] = ReconL_Con[N_x];
      //ReconR_Con[N_x] = ReconR_Con[0];
      //ReconL_Pri[0] = ReconL_Pri[N_x];
      //ReconR_Pri[N_x] = ReconR_Pri[0];
    //}
  //}
  //else if(is_RECON == 2) {//sol is primitive variables
    //VEC h(N_x);
//#pragma omp parallel for num_threads(Nthread)
    //for(u_int i = 0; i < N_x; ++i) {
      //h[i] = mesh[i+1] - mesh[i];
    //}
//#pragma omp parallel for num_threads(Nthread)
    //for(u_int i = 0; i < N_x; ++i) {//number of the cell
      //for(u_int d = 0; d < 3; ++d) {
        //if(i == 0) {
          //if(BD == 2) {
            //ENO2(h[0], h[i], h[i+1], sol[0][d], sol[i][d], sol[i+1][d], ReconR_Pri[i][d], ReconL_Pri[i+1][d]);
          //}
        //}
        //else if(i == N_x-1) {
          //if(BD == 2) {
            //ENO2(h[i-1], h[i], h[N_x-1], sol[i-1][d], sol[i][d], sol[N_x-1][d], ReconR_Pri[i][d], ReconL_Pri[i+1][d]);
          //}
        //}
        //else {
          //ENO2(h[i-1], h[i], h[i+1], sol[i-1][d], sol[i][d], sol[i+1][d], ReconR_Pri[i][d], ReconL_Pri[i+1][d]);
        //}
      //}
      //ReconL_Con[i+1] = Pri2Con(ReconL_Pri[i+1], Gamma[i]);
      //ReconR_Con[i] = Pri2Con(ReconR_Pri[i], Gamma[i]);
    //}
    //if(BD == 2) {//period
      //ReconL_Con[0] = ReconR_Con[0];
      //ReconR_Con[N_x] = ReconL_Con[N_x];
      //ReconL_Pri[0] = ReconR_Pri[0];
      //ReconR_Pri[N_x] = ReconL_Pri[N_x];
    //}
  //}

//}

//WENO
void Eulerian1D::Reconstruction(const Sol& sol, const VEC& mesh,
    Sol& ReconL_Con, Sol& ReconR_Con, Sol& ReconL_Pri, Sol& ReconR_Pri) {
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

