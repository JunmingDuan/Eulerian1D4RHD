#include "Eulerian1D.h"

bU Eulerian1D::HLL(const bU& CONL, const bU& CONR, const bU& PRIL, const bU& PRIR, const double Gammal, const double Gammar) {
  double SL, SR;//, ui, ci, vl, vr;
  double hl, csl, ul, hr, csr, ur;
  double lam1, lam2, lam3, lam4;

  //characteristic speed of left and right side
  hl = 1+PRIL[2]/PRIL[0]*Gammal/(Gammal-1);
  csl = sqrt(Gammal*PRIL[2]/PRIL[0]/hl);
  ul = PRIL[1];
  lam1 = ul-(1-ul*ul)*csl/(1-ul*csl);//negative speed
  lam2 = ul+(1-ul*ul)*csl/(1+ul*csl);//positive speed
  hr = 1+PRIR[2]/PRIR[0]*Gammar/(Gammar-1);
  csr = sqrt(Gammar*PRIR[2]/PRIR[0]/hr);
  ur = PRIR[1];
  lam3 = ur-(1-ur*ur)*csr/(1-ur*csr);//negative speed
  lam4 = ur+(1-ur*ur)*csr/(1+ur*csr);//positive speed
  SL = std::min((lam1), (lam3));//left characteristic speed
  SR = std::max((lam2), (lam4));//right characteristic speed
  if(SL >= 0) {
    return F(CONL, PRIL);
  }
  else if(SR <= 0) {
    return F(CONR, PRIR);
  }
  else {
    bU FL = F(CONL, PRIL);
    bU FR = F(CONR, PRIR);
    return (SR*FL - SL*FR + SL*SR*(CONR - CONL)) / (SR-SL);
  }
}

void Eulerian1D::cal_flux_HLL(Sol& ReconL_Con, Sol& ReconR_Con, Sol& ReconL_Pri, Sol& ReconR_Pri,
    Sol& FLUX) {
#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x+1; ++i) {
    if(i == 0) FLUX[i] = HLL(ReconL_Con[i], ReconR_Con[i], ReconL_Pri[i], ReconR_Pri[i],
        Gamma[0], Gamma[i]);
    else if(i == N_x) FLUX[i] = HLL(ReconL_Con[i], ReconR_Con[i], ReconL_Pri[i], ReconR_Pri[i],
        Gamma[i-1], Gamma[N_x-1]);
    else FLUX[i] = HLL(ReconL_Con[i], ReconR_Con[i], ReconL_Pri[i], ReconR_Pri[i],
        Gamma[i-1], Gamma[i]);
  }
}

