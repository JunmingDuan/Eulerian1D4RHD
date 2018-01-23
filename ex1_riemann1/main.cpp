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
#include "Eulerian1D.h"

bU initial(const double t, const double x, double& Gamma) {
  bU v;
  if(x < 0.5) {
    //v[0] = 10;
    //v[1] = 0;
    //v[2] = 40./3;
    //Gamma = 5./3;
    v[0] = 1;
    v[1] = 0.9;
    v[2] = 1;
    Gamma = 4./3;
  }
  else {
    //v[0] = 1;
    //v[1] = 0;
    //v[2] = 1e-6;
    //Gamma = 5./3;
    v[0] = 1;
    v[1] = 0;
    v[2] = 10;
    Gamma = 4./3;
  }

  return v;
}

int main(int argc, char* argv[])
{
  if(argc != 4)
  {
    std::cout << "Usage: " << argv[0]
      << " <Nx> <CFL> <t_end> " << std::endl;
    return 1;
  }

  double t_start = 0, t_end = atof(argv[3]);//计算时间
  double x_start = 0, x_end = 1;//计算区域
  u_int Nx = atoi(argv[1]);
  double CFL = atof(argv[2]);

  Eulerian1D Q(Nx, t_start, t_end, x_start, x_end, initial, CFL);
  std::cout << "Initialization completed ..." << std::endl;
  double t1 = omp_get_wtime();
  std::cout << "Start to solve ..." << std::endl;
  Q.Solver();
  double t2 = omp_get_wtime();

  std::ofstream outfile("sol.dat",std::ios::out);
  if(!outfile) {
    std::cout << "Open file failed!" << std::endl;
  }
  std::cout << "Print solution to " << "sol.dat" << " ..." << std::endl;
  Q.print_rupe(outfile);
  outfile.close();

  std::cout << "Time consumed: " << (t2-t1)<< std::endl;

  return 0;
}

