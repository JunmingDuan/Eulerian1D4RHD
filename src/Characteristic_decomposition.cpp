#include "Eulerian1D.h"

void Eulerian1D::cal_roe_av(const bU& PRIL, const bU& PRIR, bU& ROE_AV);
void Eulerian1D::characteristic_decomposition(const bU& CONL, const bU& CONR, bU& CHARL, bU& CHARR);
void Eulerian1D::cal_roe_av_lambda(const bU& PRIL, const bU& PRIR, double& min_lam, double& max_lam);


