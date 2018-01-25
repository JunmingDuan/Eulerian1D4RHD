#include "Eulerian1D.h"

void Eulerian1D::Solver() {
  double t_now(t_start), dt(0), alpha(0);
  Initial();
  update_cs(cs);

  int n(0);
  while(t_now < t_end) {
  //while(n < 10) {
    dt = t_step(CFL, alpha);
    dt = std::min(dt, t_end-t_now);

    //Euler_forward_LF(dt, alpha, mesh);
    //Euler_forward_LLF(dt, mesh);
    //Euler_forward_HLL(dt, mesh);
    //Euler_forward_HLLC(dt, mesh);
    //RK2_LF(Con, Pri, mesh, dt);
    //RK2_HLLC(Con, Pri, mesh, dt);
    //SSP_RK_LF(Con, Pri, mesh, dt);
    SSP_RK_HLLC(Con, Pri, mesh, dt);
    n++;

    t_now += dt;
    std::cout << "t: " << t_now << " , dt: " << dt << std::endl;
  }
  //print error of rho
  std::ofstream outfile("err.dat",std::ios::out);
  if(!outfile) {
    std::cout << "Open file failed!" << std::endl;
  }
  std::cout << "Print error to " << "err.dat" << " ..." << std::endl;
  outfile.setf(std::ios::scientific);
  outfile.precision(16);
  //outfile << "error of rho" << "\n";
  outfile << N_x << " " << cal_err(1) << " " << cal_err(2) << " " << cal_err(-1) << std::endl;
  outfile.close();

}

