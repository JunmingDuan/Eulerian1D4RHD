#include "Eulerian1D.h"

void Eulerian1D::ROE_AV_MAT(const bU& PRIL, const bU& PRIR, const double GAMMAL, const double GAMMAR,
    MAT& R, MAT& L) {
    double hl, hr, kl, kr, gl, gr, ul, ur;
    double v0, v1, v2;
    hl = 1 + GAMMAL/(GAMMAL-1)*PRIL[2]/PRIL[0];
    hr = 1 + GAMMAR/(GAMMAR-1)*PRIR[2]/PRIR[0];
    kl = sqrt(PRIL[0]*hl);
    kr = sqrt(PRIR[0]*hr);
    ul = PRIL[1];
    ur = PRIR[1];
    gl = 1./sqrt(1-ul*ul);
    gr = 1./sqrt(1-ur*ur);
    v0 = (kl*gl + kr*gr)/(kl+kr);
    v1 = (kl*gl*ul + kr*gr*ur)/(kl+kr);
    v2 = (kl*PRIL[2]/PRIL[0]/hl + kr*PRIR[2]/PRIR[0]/hr) / (kl+kr);
    double cm = 1 - GAMMAL*v2/(GAMMAL-1);
    double cp = 1 + GAMMAL*v2/(GAMMAL-1);
    double va = - v0*v0 + v1*v1;
    double s2 = 0.5*GAMMAL*v2*(1-va) - 0.5*(GAMMAL-1)*(1+va);
    double s = sqrt(s2);
    double e = -va;
    double y = sqrt((1-GAMMAL*v2)*e + s2);

    {//right characteristic matrix, R
      R[0][0] = cm;          R[0][1] = cm + s2/(GAMMAL-1); R[0][2] = cm;
      R[1][0] = v0 - s/y*v1; R[1][1] = v0;                 R[1][2] = v0 + s/y*v1;
      R[2][0] = v1 - s/y*v0; R[2][1] = v1;                 R[2][2] = v1 + s/y*v0;
    }
    {//left characteristic matrix, L
      L[0][0] = (GAMMAL-1)*e;    L[0][1] = s2*v0 - s*y*v1 - (GAMMAL-1)*e*cp*v0; L[0][2] = - s2*v1+s*y*v0 + (GAMMAL-1)*e*cp*v1;
      L[1][0] = -2*(GAMMAL-1)*e; L[1][1] = -4*s2*v0 + 2*(GAMMAL-1)*e*cp*v0;     L[1][2] = 4*s2*v1 - 2*(GAMMAL-1)*e*cp*v1;
      L[2][0] = (GAMMAL-1)*e;    L[2][1] = s2*v0 + s*y*v1 - (GAMMAL-1)*e*cp*v0; L[2][2] = - s2*v1-s*y*v0 + (GAMMAL-1)*e*cp*v1;
      L *= -0.5/e/s2;
    }
    //R[0][0] = 1; R[0][1] = 0; R[0][2] = 0;
    //R[1][0] = 0; R[1][1] = 1; R[1][2] = 0;
    //R[2][0] = 0; R[2][1] = 0; R[2][2] = 1;
    //L[0][0] = 1; L[0][1] = 0; L[0][2] = 0;
    //L[1][0] = 0; L[1][1] = 1; L[1][2] = 0;
    //L[2][0] = 0; L[2][1] = 0; L[2][2] = 1;
}

bU Eulerian1D::multiply(const bU& x, MAT& M) {
  bU y;
  y[0] = M[0][0]*x[0] + M[0][1]*x[1] + M[0][2]*x[2];
  y[1] = M[1][0]*x[0] + M[1][1]*x[1] + M[1][2]*x[2];
  y[2] = M[2][0]*x[0] + M[2][1]*x[1] + M[2][2]*x[2];
  return y;
}

