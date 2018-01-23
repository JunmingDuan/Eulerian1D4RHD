#include "Eulerian1D.h"

std::ostream& operator<<(std::ostream& os, const Eulerian1D& H) {
  os.setf(std::ios::scientific);
  os.precision(16);
  for(u_int i=0; i < H.N_x; i++) {
    os << 0.5*(H.mesh[i]+H.mesh[i+1]) << " " << H.Pri[i] << "\n";
  }
  return os;
}

void Eulerian1D::print_con(std::ostream& os) {
  os.setf(std::ios::scientific);
  os.precision(16);
  for(u_int i=0; i < N_x; i++) {
    os << 0.5*(mesh[i]+mesh[i+1]) << " " << Con[i] << "\n";
  }
}

void Eulerian1D::print_pri(std::ostream& os) {
  os.setf(std::ios::scientific);
  os.precision(16);
  for(u_int i=0; i < N_x; i++) {
    os << 0.5*(mesh[i]+mesh[i+1]) << " " << Pri[i] << "\n";
  }
}

void Eulerian1D::print_rupe(std::ostream& os) {
  os.setf(std::ios::scientific);
  os.precision(16);
  for(u_int i=0; i < N_x; i++) {
    os << 0.5*(mesh[i]+mesh[i+1]) << " "
      << Pri[i][0] << " " << Pri[i][1] << " " << Pri[i][2] << " "
      << Pri[i][2]/Pri[i][0]/(Gamma[i]-1) << "\n";
  }
}

//double Eulerian1D::fp(const bU& U, double p, const double Gamma) {
  //double gamma2 = 1 - pow(U[1]/(p+U[2]),2);
  //return U[2]+p - U[0]/sqrt(gamma2) - Gamma/(Gamma-1)*p/gamma2;
//}

/*double Eulerian1D::fpp(const bU& U, double p, const double Gamma) {*/
  //double Ep = U[2] + p;
  //double tmp = pow(U[1]/Ep, 2);
  //double gamma2 = 1 - tmp;
  //return 1 + U[0]*tmp/Ep/pow(gamma2, 1.5) - Gamma/(Gamma-1)*(gamma2 - 2*p*tmp/Ep)/pow(gamma2,2);
//}

//[>bU Eulerian1D::Con2Pri(const bU& U, const double Gamma) {<]
  ////bU prim;
  //////solve a nonlinear equation by Newton method to obtain pressure p
  ////u_int ite(0), MAXITE(10);
  ////double eps = 1e-15;
  ////double a = 0, b = (Gamma-1)*U[2];
  ////double p(0), p1(p), y(fp(U,p, Gamma));
  ////while(fabs(y) > eps && ite < MAXITE) {
    ////p1 = p - y/fpp(U, p, Gamma);
    ////y = fp(U, p1, Gamma);
    ////ite++;
    ////if(fabs(p1-p) < eps) { p = p1; break; }
    ////p = p1;
  ////}
  ////prim[2] = p;
  ////prim[1] = U[1]/(U[2]+p);
  ////prim[0] = U[0]*sqrt(1-pow(prim[1],2));

  ////return prim;
////}

//bU Eulerian1D::Con2Pri(const bU& U, const double Gamma) {
  //bU prim;
  ////solve a nonlinear equation by Newton method to obtain pressure p
  //u_int ite(0), MAXITE(20);
  //double eps = 1e-15;
  //double a = 0, b = (Gamma-1)*U[2];
  //double p(0.5*b), p1(p), y(fp(U,p, Gamma));
  //while(fabs(y) > eps && ite < MAXITE) {
    //p1 = p - y/fpp(U, p, Gamma);
    //y = fp(U, p1, Gamma);
    //ite++;
    //if(fabs(p1-p) < eps) { p = p1; break; }
    //p = p1;
  //}
  //prim[2] = p;
  //prim[1] = U[1]/(U[2]+p);
  //prim[0] = U[0]*sqrt(1-pow(prim[1],2));

  //return prim;
/*}*/

bU Eulerian1D::Pri2Con(const bU& U, const double Gamma) {
  bU Con;
  double gamma = 1./sqrt(1-U[1]*U[1]);
  double h = 1 + U[2]/U[0]*Gamma/(Gamma-1);
  Con[0] = gamma*U[0];
  Con[1] = Con[0]*h*gamma*U[1];
  Con[2] = Con[0]*h*gamma - U[2];

  return Con;
}

void Eulerian1D::update_cs(VEC& cs) {
#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x; ++i) {
    cs[i] = cal_cs(Con[i], Pri[i], Gamma[i]);
  }
}

double Eulerian1D::cal_cs
(const bU& Con, const bU& Pri, const double Gamma) {
  double h = 1+Pri[2]/Pri[0]*Gamma/(Gamma-1);
  return sqrt(Gamma*Pri[2]/Pri[0]/h);
}

double Eulerian1D::cal_max_lambda_Lag(int i) {
  double u = Pri[i][1];
  return fabs(1-u*u)*cs[i]/(1+fabs(u)*cs[i]);
}

double Eulerian1D::cal_max_lambda_Eul(int i) {
  double u = Pri[i][1];
  return fabs(1-u*u)*cs[i]/(1+fabs(u)*cs[i])+fabs(u);
}

double Eulerian1D::cal_max_lambda_Lag
(const bU& Con, const bU& Pri, const double Gamma) {
  double h = 1+Pri[2]/Pri[0]*Gamma/(Gamma-1);
  double cs = sqrt(Gamma*Pri[2]/Pri[0]/h);
  double u = Pri[1];
  return fabs(1-u*u)*cs/(1+fabs(u)*cs);
}

double Eulerian1D::cal_max_lambda_Eul
(const bU& Con, const bU& Pri, const double Gamma) {
  double h = 1+Pri[2]/Pri[0]*Gamma/(Gamma-1);
  double cs = sqrt(Gamma*Pri[2]/Pri[0]/h);
  double u = Pri[1];
  return fabs(1-u*u)*cs/(1+fabs(u)*cs)+fabs(u);
}

double Eulerian1D::cal_max_lambda_Eul() {
  double tmp, a(0);
  for(u_int i = 0; i < N_x; ++i) {
    tmp = cal_max_lambda_Eul(i);
    a = std::max(a, tmp);
  }
  return a;
}

bU Eulerian1D::F(const bU& CON, const bU& PRI) {
  bU tmp;
  tmp[0] = CON[0]*PRI[1];
  tmp[1] = CON[1]*PRI[1]+PRI[2];
  tmp[2] = CON[2]*PRI[1]+PRI[1]*PRI[2];
  return tmp;
}

double Eulerian1D::t_step(const double CFL, double& alpha) {
  double a(1), tmp_lam, hi, tmp_t;
  alpha = 0;
  for(u_int i = 0; i < N_x; ++i) {
    hi = mesh[i+1] - mesh[i];
    tmp_lam = cal_max_lambda_Eul(i);
    alpha = std::max(alpha, tmp_lam);
    tmp_t = hi/alpha;
    a = std::min(a, tmp_t);
  }
  return CFL*a;
}


