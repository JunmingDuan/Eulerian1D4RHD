#include "Eulerian1D.h"

void Eulerian1D::Euler_forward_LF(double dt, double alpha, VEC& mesh) {
  if(is_RECON == 0 || is_RECON == 1) {
    Reconstruction(Con, mesh, ReconL_Con, ReconR_Con, ReconL_Pri, ReconR_Pri);
  }
  else if(is_RECON == 2) {
    Reconstruction(Pri, mesh, ReconL_Con, ReconR_Con, ReconL_Pri, ReconR_Pri);
  }
  cal_flux_LF(ReconL_Con, ReconR_Con, ReconL_Pri, ReconR_Pri, FLUX, alpha);
#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x; ++i) {
    Con[i] = Con[i] - dt*(FLUX[i+1] - FLUX[i]) / (mesh[i+1]-mesh[i]);
  }
#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x; ++i) {
    Pri[i] = Con2Pri(Con[i], Gamma[i]);
  }
  update_cs(cs);
}

void Eulerian1D::Euler_forward_LLF(double dt, VEC& mesh) {
  if(is_RECON == 0 || is_RECON == 1) {
    Reconstruction(Con, mesh, ReconL_Con, ReconR_Con, ReconL_Pri, ReconR_Pri);
  }
  else if(is_RECON == 2) {
    Reconstruction(Pri, mesh, ReconL_Con, ReconR_Con, ReconL_Pri, ReconR_Pri);
  }
  cal_flux_LLF(ReconL_Con, ReconR_Con, ReconL_Pri, ReconR_Pri, FLUX);
#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x; ++i) {
    Con[i] = Con[i] - dt*(FLUX[i+1] - FLUX[i]) / (mesh[i+1]-mesh[i]);
  }
#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x; ++i) {
    Pri[i] = Con2Pri(Con[i], Gamma[i]);
  }
  update_cs(cs);
}

void Eulerian1D::Euler_forward_HLL(const double dt, VEC& mesh) {
  if(is_RECON == 0 || is_RECON == 1) {
    Reconstruction(Con, mesh, ReconL_Con, ReconR_Con, ReconL_Pri, ReconR_Pri);
  }
  else if(is_RECON == 2) {
    Reconstruction(Pri, mesh, ReconL_Con, ReconR_Con, ReconL_Pri, ReconR_Pri);
  }
  cal_flux_HLL(ReconL_Con, ReconR_Con, ReconL_Pri, ReconR_Pri, FLUX);
#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x; ++i) {
    Con[i] = Con[i] - dt*(FLUX[i+1] - FLUX[i]) / (mesh[i+1]-mesh[i]);
  }
#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x; ++i) {
    Pri[i] = Con2Pri(Con[i], Gamma[i]);
  }
}

void Eulerian1D::Euler_forward_HLLC(const double dt, VEC& mesh) {
  if(is_RECON == 0 || is_RECON == 1) {
    Reconstruction(Con, mesh, ReconL_Con, ReconR_Con, ReconL_Pri, ReconR_Pri);
  }
  else if(is_RECON == 2) {
    Reconstruction(Pri, mesh, ReconL_Con, ReconR_Con, ReconL_Pri, ReconR_Pri);
  }
  cal_flux_HLLC(ReconL_Con, ReconR_Con, ReconL_Pri, ReconR_Pri, FLUX);
#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x; ++i) {
    Con[i] = Con[i] - dt*(FLUX[i+1] - FLUX[i]) / (mesh[i+1]-mesh[i]);
  }
#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x; ++i) {
    Pri[i] = Con2Pri(Con[i], Gamma[i]);
  }
}

void Eulerian1D::RK2_LF(Sol& Con, Sol& Pri, VEC& mesh, const double dt) {
  Sol Con_n(Con), Pri_n(Pri);
  double alpha;
  //stage 1
  alpha = cal_max_lambda_Eul();
  if(is_RECON == 0 || is_RECON == 1) {
    Reconstruction(Con_n, mesh, ReconL_Con, ReconR_Con, ReconL_Pri, ReconR_Pri);
  }
  else if(is_RECON == 2) {
    Reconstruction(Pri_n, mesh, ReconL_Con, ReconR_Con, ReconL_Pri, ReconR_Pri);
  }
  cal_flux_LF(ReconL_Con, ReconR_Con, ReconL_Pri, ReconR_Pri, FLUX, alpha);
  //std::cout << "stage1" << std::endl;
  //std::cout << "Con\n" << std::endl;
  //for(u_int i = 0; i < N_x; ++i) {
    //std::cout << Con[i][0] << " ";
  //}
  //std::cout << "\n";
  //for(u_int i = 0; i < N_x; ++i) {
    //std::cout << Con[i][1] << " ";
  //}
  //std::cout << "\n";
  //for(u_int i = 0; i < N_x; ++i) {
    //std::cout << Con[i][2] << " ";
  //}
  //std::cout << std::endl;
  //std::cout << "ReconL_Con\n" << std::endl;
  //for(u_int i = 0; i < N_x; ++i) {
    //std::cout << ReconL_Con[i][0] << " ";
  //}
  //std::cout << "\n";
  //for(u_int i = 0; i < N_x; ++i) {
    //std::cout << ReconL_Con[i][1] << " ";
  //}
  //std::cout << "\n";
  //for(u_int i = 0; i < N_x; ++i) {
    //std::cout << ReconL_Con[i][2] << " ";
  //}
  //std::cout << std::endl;

  //std::cout << "ReconR_Con\n" << std::endl;
  //for(u_int i = 0; i < N_x; ++i) {
    //std::cout << ReconR_Con[i][0] << " ";
  //}
  //std::cout << "\n";
  //for(u_int i = 0; i < N_x; ++i) {
    //std::cout << ReconR_Con[i][1] << " ";
  //}
  //std::cout << "\n";
  //for(u_int i = 0; i < N_x; ++i) {
    //std::cout << ReconR_Con[i][2] << " ";
  //}
  //std::cout << std::endl;

#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x; ++i) {
    Con[i] = Con_n[i] - dt*(FLUX[i+1] - FLUX[i]) / (mesh[i+1]-mesh[i]);
    Pri[i] = Con2Pri(Con[i], Gamma[i]);
  }
  update_cs(cs);
  //stage 2
  alpha = cal_max_lambda_Eul();
  if(is_RECON == 0 || is_RECON == 1) {
    Reconstruction(Con, mesh, ReconL_Con, ReconR_Con, ReconL_Pri, ReconR_Pri);
  }
  else if(is_RECON == 2) {
    Reconstruction(Pri, mesh, ReconL_Con, ReconR_Con, ReconL_Pri, ReconR_Pri);
  }
  cal_flux_LF(ReconL_Con, ReconR_Con, ReconL_Pri, ReconR_Pri, FLUX, alpha);
  //std::cout << "stage2" << std::endl;
  //std::cout << "Con\n" << std::endl;
  //for(u_int i = 0; i < N_x; ++i) {
    //std::cout << Con[i][0] << " ";
  //}
  //std::cout << "\n";
  //for(u_int i = 0; i < N_x; ++i) {
    //std::cout << Con[i][1] << " ";
  //}
  //std::cout << "\n";
  //for(u_int i = 0; i < N_x; ++i) {
    //std::cout << Con[i][2] << " ";
  //}
  //std::cout << std::endl;
  //std::cout << "ReconL_Con\n" << std::endl;
  //for(u_int i = 0; i < N_x; ++i) {
    //std::cout << ReconL_Con[i][0] << " ";
  //}
  //std::cout << "\n";
  //for(u_int i = 0; i < N_x; ++i) {
    //std::cout << ReconL_Con[i][1] << " ";
  //}
  //std::cout << "\n";
  //for(u_int i = 0; i < N_x; ++i) {
    //std::cout << ReconL_Con[i][2] << " ";
  //}
  //std::cout << std::endl;

  //std::cout << "ReconR_Con\n" << std::endl;
  //for(u_int i = 0; i < N_x; ++i) {
    //std::cout << ReconR_Con[i][0] << " ";
  //}
  //std::cout << "\n";
  //for(u_int i = 0; i < N_x; ++i) {
    //std::cout << ReconR_Con[i][1] << " ";
  //}
  //std::cout << "\n";
  //for(u_int i = 0; i < N_x; ++i) {
    //std::cout << ReconR_Con[i][2] << " ";
  //}
  //std::cout << std::endl;

#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x; ++i) {
    Con[i] =  0.5*Con_n[i] + 0.5*Con[i] - 0.5*dt*(FLUX[i+1] - FLUX[i]) / (mesh[i+1]-mesh[i]);
    Pri[i] = Con2Pri(Con[i], Gamma[i]);
  }
  update_cs(cs);

}

void Eulerian1D::RK2_HLLC(Sol& Con, Sol& Pri, VEC& mesh, const double dt) {
  Sol Con_n(Con), Pri_n(Pri);
  double alpha;
  //stage 1
  if(is_RECON == 0 || is_RECON == 1) {
    Reconstruction(Con_n, mesh, ReconL_Con, ReconR_Con, ReconL_Pri, ReconR_Pri);
  }
  else if(is_RECON == 2) {
    Reconstruction(Pri_n, mesh, ReconL_Con, ReconR_Con, ReconL_Pri, ReconR_Pri);
  }
  cal_flux_HLLC(ReconL_Con, ReconR_Con, ReconL_Pri, ReconR_Pri, FLUX);
#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x; ++i) {
    Con[i] = Con_n[i] - dt*(FLUX[i+1] - FLUX[i]) / (mesh[i+1]-mesh[i]);
    Pri[i] = Con2Pri(Con[i], Gamma[i]);
  }
  update_cs(cs);
  //stage 2
  if(is_RECON == 0 || is_RECON == 1) {
    Reconstruction(Con, mesh, ReconL_Con, ReconR_Con, ReconL_Pri, ReconR_Pri);
  }
  else if(is_RECON == 2) {
    Reconstruction(Pri, mesh, ReconL_Con, ReconR_Con, ReconL_Pri, ReconR_Pri);
  }
  cal_flux_HLLC(ReconL_Con, ReconR_Con, ReconL_Pri, ReconR_Pri, FLUX);
#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x; ++i) {
    Con[i] =  0.5*Con_n[i] + 0.5*Con[i] - 0.5*dt*(FLUX[i+1] - FLUX[i]) / (mesh[i+1]-mesh[i]);
    Pri[i] = Con2Pri(Con[i], Gamma[i]);
  }
  update_cs(cs);

}

void Eulerian1D::SSP_RK_LF(Sol& Con, Sol& Pri, VEC& mesh, const double dt) {
  Sol Con_n(Con), Pri_n(Pri);
  double alpha;
  //stage 1
  alpha = cal_max_lambda_Eul();
  if(is_RECON == 0 || is_RECON == 1) {
    Reconstruction(Con_n, mesh, ReconL_Con, ReconR_Con, ReconL_Pri, ReconR_Pri);
  }
  else if(is_RECON == 2) {
    Reconstruction(Pri_n, mesh, ReconL_Con, ReconR_Con, ReconL_Pri, ReconR_Pri);
  }
  cal_flux_LF(ReconL_Con, ReconR_Con, ReconL_Pri, ReconR_Pri, FLUX, alpha);
#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x; ++i) {
    Con[i] = Con_n[i] - dt*(FLUX[i+1] - FLUX[i]) / (mesh[i+1]-mesh[i]);
    Pri[i] = Con2Pri(Con[i], Gamma[i]);
  }
  update_cs(cs);
  //stage 2
  alpha = cal_max_lambda_Eul();
  if(is_RECON == 0 || is_RECON == 1) {
    Reconstruction(Con, mesh, ReconL_Con, ReconR_Con, ReconL_Pri, ReconR_Pri);
  }
  else if(is_RECON == 2) {
    Reconstruction(Pri, mesh, ReconL_Con, ReconR_Con, ReconL_Pri, ReconR_Pri);
  }
  cal_flux_LF(ReconL_Con, ReconR_Con, ReconL_Pri, ReconR_Pri, FLUX, alpha);
#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x; ++i) {
    Con[i] = 0.75*Con_n[i] + 0.25*Con[i] - 0.25*dt*(FLUX[i+1] - FLUX[i]) / (mesh[i+1]-mesh[i]);
    Pri[i] = Con2Pri(Con[i], Gamma[i]);
  }
  update_cs(cs);
  //stage 3
  alpha = cal_max_lambda_Eul();
  if(is_RECON == 0 || is_RECON == 1) {
    Reconstruction(Con, mesh, ReconL_Con, ReconR_Con, ReconL_Pri, ReconR_Pri);
  }
  else if(is_RECON == 2) {
    Reconstruction(Pri, mesh, ReconL_Con, ReconR_Con, ReconL_Pri, ReconR_Pri);
  }
  cal_flux_LF(ReconL_Con, ReconR_Con, ReconL_Pri, ReconR_Pri, FLUX, alpha);
#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x; ++i) {
    Con[i] = 1./3*Con_n[i] + 2./3*Con[i] - 2./3*dt*(FLUX[i+1] - FLUX[i]) / (mesh[i+1]-mesh[i]);
    Pri[i] = Con2Pri(Con[i], Gamma[i]);
  }
  update_cs(cs);

}

void Eulerian1D::SSP_RK_HLLC(Sol& Con, Sol& Pri, VEC& mesh, const double dt) {
  Sol Con_n(Con), Pri_n(Pri);
  //stage 1
  if(is_RECON == 0 || is_RECON == 1) {
    Reconstruction(Con_n, mesh, ReconL_Con, ReconR_Con, ReconL_Pri, ReconR_Pri);
  }
  else if(is_RECON == 2) {
    Reconstruction(Pri_n, mesh, ReconL_Con, ReconR_Con, ReconL_Pri, ReconR_Pri);
  }
  cal_flux_HLLC(ReconL_Con, ReconR_Con, ReconL_Pri, ReconR_Pri, FLUX);
#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x; ++i) {
    Con[i] = Con_n[i] - dt*(FLUX[i+1] - FLUX[i]) / (mesh[i+1]-mesh[i]);
    Pri[i] = Con2Pri(Con[i], Gamma[i]);
  }
  update_cs(cs);
  //stage 2
  if(is_RECON == 0 || is_RECON == 1) {
    Reconstruction(Con, mesh, ReconL_Con, ReconR_Con, ReconL_Pri, ReconR_Pri);
  }
  else if(is_RECON == 2) {
    Reconstruction(Pri, mesh, ReconL_Con, ReconR_Con, ReconL_Pri, ReconR_Pri);
  }
  cal_flux_HLLC(ReconL_Con, ReconR_Con, ReconL_Pri, ReconR_Pri, FLUX);
#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x; ++i) {
    Con[i] = 0.75*Con_n[i] + 0.25*Con[i] - 0.25*dt*(FLUX[i+1] - FLUX[i]) / (mesh[i+1]-mesh[i]);
    Pri[i] = Con2Pri(Con[i], Gamma[i]);
  }
  update_cs(cs);
  //stage 3
  if(is_RECON == 0 || is_RECON == 1) {
    Reconstruction(Con, mesh, ReconL_Con, ReconR_Con, ReconL_Pri, ReconR_Pri);
  }
  else if(is_RECON == 2) {
    Reconstruction(Pri, mesh, ReconL_Con, ReconR_Con, ReconL_Pri, ReconR_Pri);
  }
  cal_flux_HLLC(ReconL_Con, ReconR_Con, ReconL_Pri, ReconR_Pri, FLUX);
#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x; ++i) {
    Con[i] = 1./3*Con_n[i] + 2./3*Con[i] - 2./3*dt*(FLUX[i+1] - FLUX[i]) / (mesh[i+1]-mesh[i]);
    Pri[i] = Con2Pri(Con[i], Gamma[i]);
  }
  update_cs(cs);
}
  //std::cout << "stage1" << std::endl;
  //std::cout << "Con\n" << std::endl;
  //for(u_int i = 0; i < N_x; ++i) {
    //std::cout << Con[i][0] << " ";
  //}
  //std::cout << "\n";
  //for(u_int i = 0; i < N_x; ++i) {
    //std::cout << Con[i][1] << " ";
  //}
  //std::cout << "\n";
  //for(u_int i = 0; i < N_x; ++i) {
    //std::cout << Con[i][2] << " ";
  //}
  //std::cout << std::endl;
  //std::cout << "ReconL_Con\n" << std::endl;
  //for(u_int i = 0; i < N_x+1; ++i) {
    //std::cout << ReconL_Con[i][0] << " ";
  //}
  //std::cout << "\n";
  //for(u_int i = 0; i < N_x+1; ++i) {
    //std::cout << ReconL_Con[i][1] << " ";
  //}
  //std::cout << "\n";
  //for(u_int i = 0; i < N_x+1; ++i) {
    //std::cout << ReconL_Con[i][2] << " ";
  //}
  //std::cout << std::endl;

  //std::cout << "ReconR_Con\n" << std::endl;
  //for(u_int i = 0; i < N_x+1; ++i) {
    //std::cout << ReconR_Con[i][0] << " ";
  //}
  //std::cout << "\n";
  //for(u_int i = 0; i < N_x+1; ++i) {
    //std::cout << ReconR_Con[i][1] << " ";
  //}
  //std::cout << "\n";
  //for(u_int i = 0; i < N_x+1; ++i) {
    //std::cout << ReconR_Con[i][2] << " ";
  //}
  //std::cout << std::endl;


