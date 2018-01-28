#include "Eulerian1D.h"
#include "WENO_nonuniform.h"
#include "ENO.h"

//ENO
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
      ReconL_Con[0] = ReconR_Con[N_x];
      ReconR_Con[N_x] = ReconL_Con[0];
      ReconL_Pri[0] = ReconR_Pri[N_x];
      ReconR_Pri[N_x] = ReconL_Pri[0];
    }
  }

}

////WENO
//void Eulerian1D::Reconstruction(const Sol& sol, const VEC& mesh,
    //Sol& ReconL_Con, Sol& ReconR_Con, Sol& ReconL_Pri, Sol& ReconR_Pri) {
  //if(is_RECON == 0) {//sol,gl,gr are conservative variables
//#pragma omp parallel for num_threads(Nthread)
    //for(u_int i = 0; i < N_x+1; ++i) {
      //if(i == 0) {
        //if(BD == 1) {
          //ReconL_Con[i] = sol[0];
        //}
        //else if(BD == 2) {
          //ReconL_Con[i] = sol[N_x-1];
        //}
        //ReconR_Con[i] = sol[i];
        //ReconL_Pri[i] = Con2Pri(ReconL_Con[i], Gamma[N_x-1]);
        //ReconR_Pri[i] = Con2Pri(ReconR_Con[i], Gamma[i]);
      //}
      //else if(i == N_x) {
        //if(BD == 1) {
          //ReconR_Con[i] = sol[N_x-1];
        //}
        //else if(BD == 2) {
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
          //if(BD == 1) {
            //WENO3(h[0], h[i], h[i+1], sol[0][d], sol[i][d], sol[i+1][d], ReconR_Con[i][d], ReconL_Con[i+1][d]);
          //}
          //else if(BD == 2) {
            //WENO3(h[N_x-1], h[i], h[i+1], sol[N_x-1][d], sol[i][d], sol[i+1][d], ReconR_Con[i][d], ReconL_Con[i+1][d]);
          //}
        //}
        //else if(i == N_x-1) {
/*bU Lagranian1D::HLLC(const bU& CONL, const bU& CONR, const bU& PRIL, const bU& PRIR, const double Gammal, const double Gammar, double& SM) {*/
  //double SL, SR, ui, ci, vl, vr;
  //double srhol, srhor;
  //double hl, csl, ul, hr, csr, ur;
  //double lam1, lam2, lam3, lam4;

  ////characteristic speed of left and right side
  //hl = 1+PRIL[2]/PRIL[0]*Gammal/(Gammal-1);
  //csl = sqrt(Gammal*PRIL[2]/PRIL[0]/hl);
  //ul = PRIL[1];
  //lam1 = ul-(1-ul*ul)*csl/(1-ul*csl);//negative speed
  //lam2 = ul+(1-ul*ul)*csl/(1+ul*csl);//positive speed
  //hr = 1+PRIR[2]/PRIR[0]*Gammar/(Gammar-1);
  //csr = sqrt(Gammar*PRIR[2]/PRIR[0]/hr);
  //ur = PRIR[1];
  //lam3 = ur-(1-ur*ur)*csr/(1-ur*csr);//negative speed
  //lam4 = ur+(1-ur*ur)*csr/(1+ur*csr);//positive speed
  ////roe_average of ul, ur
  //srhol = sqrt(PRIL[0]);
  //srhor = sqrt(PRIR[0]);
  //ui = (srhol*ul + srhor*ur)/(srhol+srhor);
  ////left and right speed limit
  //vl = ul - ui;
  //vr = ur - ui;
  //SL = lam1 - ui;//left characteristic speed
  //SR = lam4 - ui;//right characteristic speed

  ////S1 = std::min(vl-cl, -ci);
  ////S2 = std::max(vr+cr, ci);

  //if(SL > 0) {
    //SM = ui;
    //return F(CONL, PRIL);
  //}
  //else if(SR < 0) {
    //SM = ui;
    //return F(CONR, PRIR);
  //}
  //else
  //{
    ////cal S^* and p^*, S^* is the smaller root of a quadric equation:
    ////F_E^{hll} x^2 - (E^{hll}+F^{hll}_m) x + m^{hll} = 0
    //double coe1 = SR*CONL[1] - SL*CONR[1] + SL*SR*(CONR[2] - CONL[2]);
    //double coe2 = ( SR*(CONL[1]*PRIL[1]+PRIL[2]) - SL*(CONR[1]*PRIR[1]+PRIR[2]) + SL*SR*(CONR[1] - CONL[1])
        //+ SR*CONR[2] - SL*CONL[2] + CONL[1] - CONR[1] );
    //double coe3 = SR*CONR[1] - SL*CONL[1] + CONL[1]*PRIL[1]+PRIL[2] - CONR[1]*PRIR[1]-PRIR[2];
    ////std::cout << "SL,SR: " << SL << " " << SR << std::endl;
    ////std::cout << "COE: " << coe1 << " " << coe2 << " " << coe3 << std::endl;
    //double PM;
    //if(fabs(coe1) < 1e-10) SM = coe3/coe2;
          //if(BD == 1) {
            //WENO3(h[i-1], h[i], h[N_x-1], sol[i-1][d], sol[i][d], sol[N_x-1][d], ReconR_Con[i][d], ReconL_Con[i+1][d]);
          //}
          //else if(BD == 2) {
            //WENO3(h[i-1], h[i], h[0], sol[i-1][d], sol[i][d], sol[0][d], ReconR_Con[i][d], ReconL_Con[i+1][d]);
          //}
        //}
        //else {
          //WENO3(h[i-1], h[i], h[i+1], sol[i-1][d], sol[i][d], sol[i+1][d], ReconR_Con[i][d], ReconL_Con[i+1][d]);
        //}
      //}
      //ReconL_Pri[i+1] = Con2Pri(ReconL_Con[i+1], Gamma[i]);
      //ReconR_Pri[i] = Con2Pri(ReconR_Con[i], Gamma[i]);
    //}
    //if(BD == 1) {//outflow
      //ReconL_Con[0] = ReconR_Con[0];
      //ReconR_Con[N_x] = ReconL_Con[N_x];
      //ReconL_Pri[0] = ReconR_Pri[0];
      //ReconR_Pri[N_x] = ReconL_Pri[N_x];
    //}
    //else if(BD == 2) {//period
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
          //if(BD == 1) {
            //WENO3(h[0], h[i], h[i+1], sol[0][d], sol[i][d], sol[i+1][d], ReconR_Pri[i][d], ReconL_Pri[i+1][d]);
          //}
          //else if(BD == 2) {
            //WENO3(h[N_x-1], h[i], h[i+1], sol[N_x-1][d], sol[i][d], sol[i+1][d], ReconR_Pri[i][d], ReconL_Pri[i+1][d]);
          //}
        //}
        //else if(i == N_x-1) {
          //if(BD == 1) {
            //WENO3(h[i-1], h[i], h[N_x-1], sol[i-1][d], sol[i][d], sol[N_x-1][d], ReconR_Pri[i][d], ReconL_Pri[i+1][d]);
          //}
          //else if(BD == 2) {
            //WENO3(h[i-1], h[i], h[0], sol[i-1][d], sol[i][d], sol[0][d], ReconR_Pri[i][d], ReconL_Pri[i+1][d]);
          //}
        //}
        //else {
          //WENO3(h[i-1], h[i], h[i+1], sol[i-1][d], sol[i][d], sol[i+1][d], ReconR_Pri[i][d], ReconL_Pri[i+1][d]);
        //}
      //}
      //ReconL_Con[i+1] = Pri2Con(ReconL_Pri[i+1], Gamma[i]);
      //ReconR_Con[i] = Pri2Con(ReconR_Pri[i], Gamma[i]);
    //}
    //if(BD == 1) {//outflow
      //ReconL_Con[0] = ReconR_Con[0];
      //ReconR_Con[N_x] = ReconL_Con[N_x];
      //ReconL_Pri[0] = ReconR_Pri[0];
      //ReconR_Pri[N_x] = ReconL_Pri[N_x];
    //}
    //else if(BD == 2) {//period
      //ReconL_Con[0] = ReconL_Con[N_x];
      //ReconR_Con[N_x] = ReconR_Con[0];
      //ReconL_Pri[0] = ReconL_Pri[N_x];
      //ReconR_Pri[N_x] = ReconR_Pri[0];
    //}
  //}

//}

