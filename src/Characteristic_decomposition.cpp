#include "Eulerian1D.h"

void Eulerian1D::cal_roe_av(const bU& PRIL, const bU& PRIR, const double GammaL, const double GammaR, bU& ROE_AV) {
  double hL = 1 + PRIL[2]/PRIL[0]*GammaL/(GammaL-1);
  double hR = 1 + PRIR[2]/PRIR[0]*GammaR/(GammaR-1);
  double kL = PRIL[0]*hL;
  double kR = PRIR[0]*hR;
  double gammaL = sqrt(1-PRIL[1]*PRIL[1]);
  double gammaR = sqrt(1-PRIR[1]*PRIR[1]);
  double wL, wR;
  wL = kL*gammaL;
  wR = kR*gammaR;
  ROE_AV[0] = (wL + wR)/(kL + kR);
}
void Eulerian1D::characteristic_decomposition(const bU& CONL, const bU& CONR, bU& CHARL, bU& CHARR) {
}
void Eulerian1D::cal_roe_av_lambda(const bU& PRIL, const bU& PRIR, double& min_lam, double& max_lam) {
}

