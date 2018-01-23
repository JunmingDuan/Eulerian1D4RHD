#include "gsl/gsl_roots.h"
#include "Eulerian1D.h"

double Eulerian1D::fp(double x, void *params)
{
  struct Con_params *p = (struct Con_params *) params;

  double D = p->D;
  double m = p->m;
  double E = p->E;
  double Gamma = p->Gamma;
  double gamma2 = 1 - pow(m/(x+E),2);
  return E+x - D/sqrt(gamma2) - Gamma/(Gamma-1)*x/gamma2;
}

double Eulerian1D::fpp(double x, void *params)
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

void Eulerian1D::fp_fdf(double x, void *params, double *y, double *dy)
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

bU Eulerian1D::Con2Pri(const bU& U, const double Gamma) {
  int status;
  int iter = 0, max_iter = 20;
  const gsl_root_fdfsolver_type *T;
  gsl_root_fdfsolver *s;
  double x0, x = 0.5*(Gamma-1)*U[2];
  gsl_function_fdf FDF;
  struct Con_params params = {U[0], U[1], U[2], Gamma};

  FDF.f = &fp;
  FDF.df = &fpp;
  FDF.fdf = &fp_fdf;
  FDF.params = &params;

  T = gsl_root_fdfsolver_newton;
  s = gsl_root_fdfsolver_alloc (T);
  gsl_root_fdfsolver_set (s, &FDF, x);

  printf ("using %s method\n",
      gsl_root_fdfsolver_name (s));

  printf ("%-5s %10s %10s %10s\n",
      "iter", "root", "err", "err(est)");
  do
  {
    iter++;
    status = gsl_root_fdfsolver_iterate (s);
    x0 = x;
    x = gsl_root_fdfsolver_root (s);
    status = gsl_root_test_delta (x, x0, 0, 1e-3);

    if (status == GSL_SUCCESS)
      printf ("Converged:\n");

    printf ("%5d %10.7f %10.7f\n",
        iter, x,  x - x0);
  }
  while (status == GSL_CONTINUE && iter < max_iter);

  gsl_root_fdfsolver_free (s);

  abort();
  return x;
}


