#include "ENO.h"

template <class T>
void ENO2(double h1, double h2, double h3, const T& u1, const T& u2, const T& u3, const std::vector<double>& point,
    T& ul, T& ur) {
  double eps(1e-13);
  //rho
  //poly
  double k1 = 2.*(u2[0] - u1[0])/(h1 + h2);
  double k2 = 2.*(u3[0] - u2[0])/(h2 + h3);
  //limiter
  double theta(1.0), krho;
  if(fabs(k1) > fabs(k2)) { krho = k2; }
  else { krho = k1; }

  for(u_int i = 0; i < point.size(); ++i) {
    double up = krho*(h2*point[i] - h2/2.) + u2[0];//u[0] at the point
    if(fabs(u2[0] - up) < 0.1*eps) continue;
    else if(fabs((u2[0]-eps)/(u2[0]-up)) < theta)
      theta = fabs((u2[0]-eps)/(u2[0]-up));
  }
  ul[0] = theta*krho*(h2*point.front() - h2/2.) + u2[0];
  ur[0] = theta*krho*(h2*point.back() - h2/2.) + u2[0];

  //u
  //poly
  k1 = 2.*(u2[1] - u1[1])/(h1 + h2);
  k2 = 2.*(u3[1] - u2[1])/(h2 + h3);
  //limiter
  double ku;
  if(fabs(k1) > fabs(k2)) { ku = k2; }
  else { ku = k1; }

  ul[1] = ku*(h2*point.front() - h2/2.) + u2[1];
  ur[1] = ku*(h2*point.back() - h2/2.) + u2[1];

  //E
  //poly
  k1 = 2.*(u2[2] - u1[2])/(h1 + h2);
  k2 = 2.*(u3[2] - u2[2])/(h2 + h3);
  //limiter
  theta = 1.0;
  double kE;
  if(fabs(k1) > fabs(k2)) { kE = k2; }
  else { kE = k1; }

  std::vector<double> e(point.size());
  int flag(1);
  double ec = u2[2] - 0.5*pow(u2[1], 2);
  for(u_int i = 0; i < point.size(); ++i) {
    e[i] = kE*(h2*point[i] - h2/2.) + u2[2] - 0.5*pow(ku*(h2*point[i] - h2/2.) + u2[1] ,2);//e[2] at the point
  }
  for(u_int i = 0; i < point.size(); ++i) {
    if(e[i] < eps) { flag = 0; break; }
  }
  if(flag == 1) {
    theta = 1.0;
  }
  else {
    for(u_int i = 0; i < point.size(); ++i) {
      if(ec/(ec-e[i]) < theta)
        theta = ec/(ec-e[i]);
    }
  }

  ul[2] = theta*kE*(h2*point.front() - h2/2.) + u2[2];
  ur[2] = theta*kE*(h2*point.back() - h2/2.) + u2[2];

}

double ENO2_CELL_L(double h1, double h2, double h3, const double u1, const double u2, const double u3) {
  double k1 = 2.*(u2 - u1)/(h1 + h2);
  double k2 = 2.*(u3 - u2)/(h2 + h3);
  double krho;
  krho = fabs(k1) > fabs(k2) ? k2 : k1;

  //return u2;
  return -krho*(h2/2.) + u2;
}

double ENO2_CELL_R(double h1, double h2, double h3, const double u1, const double u2, const double u3) {
  double k1 = 2.*(u2 - u1)/(h1 + h2);
  double k2 = 2.*(u3 - u2)/(h2 + h3);
  double krho;
  krho = fabs(k1) > fabs(k2) ? k2 : k1;

  //return u2;
  return krho*(h2/2.) + u2;
}

void ENO2(double h1, double h2, double h3, const double u1, const double u2, const double u3,
    double& ul, double& ur) {
  //poly
  double k1 = 2.*(u2 - u1)/(h1 + h2);
  double k2 = 2.*(u3 - u2)/(h2 + h3);
  //minmod limiter
  double krho;
  krho = fabs(k1) > fabs(k2) ? k2 : k1;

  ul = -krho*(h2/2.) + u2;
  ur = krho*(h2/2.) + u2;
}


