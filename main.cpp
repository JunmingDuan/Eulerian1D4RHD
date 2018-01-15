/**
 * @file main.cpp
 * @brief test cpp for 1D Lagranian scheme for RHD
 * @author Duan Junming, duanjm@pku.edu.cn
 * @version 1.0
 * @date 2018-01-15
 */
#include <iostream>
#include <cstdlib>
#include <time.h>
#include <string>
#include <sstream>
#include <iomanip>
#include <fstream>
#include "vvector.h"
#include "mvector.h"
#include "Lag1D.h.h"

bU initial(const double t, const double x, double& gamma) {
  bU v;
  //if(x < 0.1) {
    //double rhol = 1;
    //double ul = 1;
    //double pl = 1e3;
    //gamma = 1.4;
    //v[0] = rhol;
    //v[1] = ul;
    //v[2] = pl/(gamma-1.0) + 0.5*ul*ul*rhol;
  //}
  //else if(x < 0.9) {
    //double rhom = 1;
    //double um = 1;
    //double pm = 1e-2;
    //gamma = 1.4;
    //v[0] = rhom;
    //v[1] = um;
    //v[2] = pm/(gamma-1.0) + 0.5*um*um*rhom;
  //}
  //else {
    //double rhor = 1;
    //double ur = 1;
    //double pr = 1e2;
    //gamma = 1.4;
    //v[0] = rhor;
    //v[1] = ur;
    //v[2] = pr/(gamma-1.0) + 0.5*ur*ur*rhor;
  //}
  double p = 1e-13;
  gamma = 5./3;
  v[0] = 1;
  v[1] = -1;
  v[2] = p/(gamma-1.0)+0.5;

  return v;
}

int main(int argc, char* argv[])
{
  if(argc < 4)
  {
    std::cout << "Usage: " << argv[0]
      << " <Nx> <CFL> <t_end> " << std::endl;
    return 1;
  }

  double t_start = 0,  t_end = atof(argv[3]);//计算时间
  double x_start = 0, x_end = 1;//计算区域
  u_int Nx = atoi(argv[1]);
  double CFL = atof(argv[2]);

  SCL1D<bU> Q(Nx, t_start, t_end, x_start, x_end,
      initial, CFL);
  clock_t t1, t2;
  t1 = clock();
  Q.Solve();
  t2 = clock();

  std::ofstream outfile("sol.dat",std::ios::out);
  if(!outfile) {
    std::cout << "Open file failed!" << std::endl;
  }

  outfile << Q << std::endl;
  std::cout << "time: " << (double)(t2-t1)/CLOCKS_PER_SEC << std::endl;

  return 0;
}

