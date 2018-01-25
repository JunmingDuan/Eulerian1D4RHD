#include "Eulerian1D.h"

double func_p(double x, void *params)
{
  struct Con_params *p = (struct Con_params *) params;

  double D = p->D;
  double m = p->m;
  double E = p->E;
  double Gamma = p->Gamma;
  double gamma2 = 1 - pow(m/(x+E),2);
  return E+x - D/sqrt(gamma2) - Gamma/(Gamma-1)*x/gamma2;
}

double func_dp(double x, void *params)
{
  struct Con_params *p = (struct Con_params *) params;
  double D = p->D;
  double m = p->m;
  double E = p->E;
  double Gamma = p->Gamma;
  double Ep = E + x;
  double tmp = pow(m/Ep, 2);
  double gamma2 = 1 - tmp;
  return 1 + D*tmp/Ep/pow(gamma2, 1.5) - Gamma/(Gamma-1)*(gamma2 - 2*x*tmp/Ep)/pow(gamma2,2);
}

void func_pdp(double x, void *params, double *y, double *dy)
{
  struct Con_params *p = (struct Con_params *) params;
  double D = p->D;
  double m = p->m;
  double E = p->E;
  double Gamma = p->Gamma;
  double Ep = E + x;
  double tmp = pow(m/Ep, 2);
  double gamma2 = 1 - tmp;
  *y = E+x - D/sqrt(gamma2) - Gamma/(Gamma-1)*x/gamma2;
  *dy = 1 + D*tmp/Ep/pow(gamma2, 1.5) - Gamma/(Gamma-1)*(gamma2 - 2*x*tmp/Ep)/pow(gamma2,2);
}

double Eulerian1D::fp(const bU& U, double p, const double Gamma) {
  double gamma2 = 1 - pow(U[1]/(p+U[2]),2);
  return U[2]+p - U[0]/sqrt(gamma2) - Gamma/(Gamma-1)*p/gamma2;
}

double Eulerian1D::fpp(const bU& U, double p, const double Gamma) {
  double Ep = U[2] + p;
  double tmp = pow(U[1]/Ep, 2);
  double gamma2 = 1 - tmp;
  return 1 + U[0]*tmp/Ep/pow(gamma2, 1.5) - Gamma/(Gamma-1)*(gamma2 - 2*p*tmp/Ep)/pow(gamma2,2);
}

//double Eulerian1D::fp(const bU& U, double p, const double Gamma) {
  //double Ep = U[2]+p;
  //double u = U[1]/Ep;
  //double gamma2 = 1-u*u;
  //return Ep*gamma2 - U[0]*sqrt(gamma2) - Gamma/(Gamma-1)*p;
//}

//double Eulerian1D::fpp(const bU& U, double p, const double Gamma) {
  //double Ep = U[2]+p;
  //double u = U[1]/Ep;
  //double gamma2 = 1-u*u;
  //return 1 - 3*u*u - Gamma/(Gamma-1)*p - U[0]*pow(u,3)/U[1]/sqrt(gamma2);
//}

//solve a nonlinear equation by Newton method to obtain pressure p
bU Eulerian1D::Con2Pri(const bU& U, const double Gamma) {
  bU prim;
  u_int ite(0), MAXITE(10);
  double eps = 1e-15;
  double a = 0, b = (Gamma-1)*U[2];
  double p(0.5*b), p1(p), y(fp(U,p, Gamma));
  //std::cout.setf(std::ios::scientific);
  //std::cout.precision(16);
  //std::cout << "E: " << U[2] << std::endl;
  //std::cout << "b: " << b << std::endl;
  while(fabs(y) > eps && ite < MAXITE) {
    p1 = p - y/fpp(U, p, Gamma);
    if(fabs(p1-p) < eps) { p = p1; break; }
    //std::cout << "ite: " << ite << " ; p1: " << p1 << std::endl;
    //while(p1 < a) {
      //p1 = 0.5*(p1+b);
      //std::cout << "<a: " << p1 << std::endl;
    //}
    //while(p1 > b) {
      //p1 = 0.5*(p1+a);
      //std::cout << ">b: " << p1 << std::endl;
    //}
    if(p1 < a) p1 = a;
    if(p1 > b) p1 = b;
    y = fp(U, p1, Gamma);
    //std::cout << "ite: " << ite << " ; corrected p1: " << p1 << std::endl;
    //std::cout << "y: " << y << std::endl;
    p = p1;
    ite++;
  }
  if(p < 0 || p != p) {
    std::cout << U[0] << " " << U[1] << " " << U[2] << std::endl;
    std::cout << "ite: " << ite << " ; p: " << p << " ; y: " << y << std::endl;
    std::cout << "derivative: " << fpp(U, 0.5*b, Gamma) << std::endl;
    //std::cout << "p1: " << 0.5*b - fp(U,0.5*b,Gamma)/fpp(U, 0.5*b, Gamma) << std::endl;
    std::cout << "f(p1): " << fp(U,p, Gamma) << std::endl;
    std::cout << U[1]/(U[2]+p) << std::endl;
    abort();
  }

  prim[2] = p;
  prim[1] = U[1]/(U[2]+p);
  prim[0] = U[0]*sqrt(1-pow(prim[1],2));
  //std::cout << "ite: " << ite << " ; final p: " << p << std::endl;
  //std::cout << "rho: " << prim[0] << std::endl;
  //abort();

  return prim;
}

//bU Eulerian1D::Con2Pri(const bU& U, const double Gamma) {
  //int status1, status2;
  //int iter = 0;
  //double x0, x = 0.5*(Gamma-1)*U[2];
  //struct Con_params params = {U[0], U[1], U[2], Gamma};
  //double y;

  //FDF.params = &params;
  //gsl_root_fdfsolver_set (s, &FDF, x);

  //std::cout << U[0] << " " << U[1] << " " << U[2] << std::endl;
  //do
  //{
    //iter++;
    //status1 = gsl_root_fdfsolver_iterate (s);
    //x0 = x;
    //x = gsl_root_fdfsolver_root (s);
    //status1 = gsl_root_test_delta (x, x0, 1e-10, 1e-10);
    //y = GSL_FN_FDF_EVAL_F(&FDF, x);
    //status2 = gsl_root_test_residual(y, 1e-10);
    //std::cout << "iter: " << iter  << " ; x: " << x << std::endl;
  //}
  //while (status1 == GSL_CONTINUE && status2 == GSL_CONTINUE && iter < 10);
  //if(x < 0) {
    //std::cout << x << std::endl;
    //abort();
  //}
  //bU U1;
  //U1[2] = x;
  //U1[1] = U[1]/(U[2]+x);
  //U1[0] = U[0]*sqrt(1-pow(U1[1],2));

  //return U1;
//}


