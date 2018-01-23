#include "Eulerian1D.h"

void Eulerian1D::Euler_forward_LF(double dt, double alpha, VEC& mesh) {
  Reconstruction(Con, mesh, ReconL_Con, ReconR_Con, ReconL_Pri, ReconR_Pri);
  cal_flux_LF(ReconL_Con, ReconR_Con, ReconL_Pri, ReconR_Pri, FLUX, alpha);
  cal_us_roeav(ReconL_Pri, ReconR_Pri, us);
#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x; ++i) {
    Con[i] *= (mesh[i+1]-mesh[i]);
  }
  move_mesh(mesh, us, dt, mesh);
#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x; ++i) {
    Con[i] = (Con[i] - dt*(FLUX[i+1] - FLUX[i])) / (mesh[i+1]-mesh[i]);
  }
#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x; ++i) {
    Pri[i] = Con2Pri(Con[i], Gamma[i]);
  }
  update_cs(cs);
}

void Eulerian1D::Euler_forward_LLF(double dt, VEC& mesh) {
  Reconstruction(Con, mesh, ReconL_Con, ReconR_Con, ReconL_Pri, ReconR_Pri);
  cal_flux_LLF(ReconL_Con, ReconR_Con, ReconL_Pri, ReconR_Pri, FLUX);
  cal_us_roeav(ReconL_Pri, ReconR_Pri, us);
#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x; ++i) {
    Con[i] *= (mesh[i+1]-mesh[i]);
  }
  move_mesh(mesh, us, dt, mesh);
#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x; ++i) {
    Con[i] = (Con[i] - dt*(FLUX[i+1] - FLUX[i])) / (mesh[i+1]-mesh[i]);
  }
#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x; ++i) {
    Pri[i] = Con2Pri(Con[i], Gamma[i]);
  }
  update_cs(cs);
}

void Eulerian1D::RK2_LF(Sol& Con, Sol& Pri, VEC& mesh, const double dt) {
  VEC mesh_n(mesh), mesh1(N_x+1);
  Sol Con_n(Con), Pri_n(Pri);
  double alpha;
  //stage 1
  alpha = cal_max_lambda_Eul();
  Reconstruction(Con_n, mesh_n, ReconL_Con, ReconR_Con, ReconL_Pri, ReconR_Pri);
  //Reconstruction(Pri, mesh, ReconL_Con, ReconR_Con, ReconL_Pri, ReconR_Pri);
  cal_flux_LF(ReconL_Con, ReconR_Con, ReconL_Pri, ReconR_Pri, FLUX, alpha);
  cal_us_roeav(ReconL_Pri, ReconR_Pri, us);
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
  for(u_int i = 0; i < N_x+1; ++i) {
    mesh[i] = mesh_n[i] + dt * us[i];
  }
#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x; ++i) {
    Con[i] = (Con_n[i]*(mesh_n[i+1]-mesh_n[i]) - dt*(FLUX[i+1] - FLUX[i])) / (mesh[i+1]-mesh[i]);
    Pri[i] = Con2Pri(Con[i], Gamma[i]);
  }
  update_cs(cs);
  //stage 2
  alpha = cal_max_lambda_Eul();
  Reconstruction(Con, mesh, ReconL_Con, ReconR_Con, ReconL_Pri, ReconR_Pri);
  //Reconstruction(Pri, mesh, ReconL_Con, ReconR_Con, ReconL_Pri, ReconR_Pri);
  cal_flux_LF(ReconL_Con, ReconR_Con, ReconL_Pri, ReconR_Pri, FLUX, alpha);
  cal_us_roeav(ReconL_Pri, ReconR_Pri, us);
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
  for(u_int i = 0; i < N_x+1; ++i) {
    mesh1[i] = 0.5*mesh_n[i] + 0.5*(mesh[i] + dt * us[i]);
  }
#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x; ++i) {
    Con[i] = ( 0.5*Con_n[i]*(mesh_n[i+1]-mesh_n[i])
        + 0.5*(Con[i]*(mesh[i+1]-mesh[i]) - dt*(FLUX[i+1] - FLUX[i])) )/ (mesh1[i+1]-mesh1[i]);
    Pri[i] = Con2Pri(Con[i], Gamma[i]);
  }
  mesh = mesh1;
  update_cs(cs);

}


void Eulerian1D::Euler_forward_HLLC(const double dt, VEC& mesh) {
  Reconstruction(Con, mesh, ReconL_Con, ReconR_Con, ReconL_Pri, ReconR_Pri);
  cal_flux_HLLC(ReconL_Con, ReconR_Con, ReconL_Pri, ReconR_Pri, FLUX, us);
#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x; ++i) {
    Con[i] *= (mesh[i+1]-mesh[i]);
  }
  move_mesh(mesh, us, dt, mesh);
#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x; ++i) {
    Con[i] = (Con[i] - dt*(FLUX[i+1] - FLUX[i])) / (mesh[i+1]-mesh[i]);
  }
#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x; ++i) {
    Pri[i] = Con2Pri(Con[i], Gamma[i]);
  }
}

void Eulerian1D::SSP_RK_LF(Sol& Con, Sol& Pri, VEC& mesh, const double dt, double alpha) {
  VEC mesh_n(mesh), mesh1(N_x+1);
  Sol Con_n(Con), Pri_n(Pri);
  //stage 1
  Reconstruction(Con, mesh, ReconL_Con, ReconR_Con, ReconL_Pri, ReconR_Pri);
  cal_flux_LF(ReconL_Con, ReconR_Con, ReconL_Pri, ReconR_Pri, FLUX, alpha);
  cal_us_roeav(ReconL_Pri, ReconR_Pri, us);
  std::cout << "stage1" << std::endl;
  std::cout << "Con\n" << Con << std::endl;
  std::cout << "Pri\n" << Pri << std::endl;
  std::cout << "ReconL_Con\n" << ReconL_Con << std::endl;
  std::cout << "ReconR_Con\n" << ReconR_Con << std::endl;
  std::cout << "ReconL_Pri\n" << ReconL_Pri << std::endl;
  std::cout << "ReconR_Pri\n" << ReconR_Pri << std::endl;
  std::cout << "FLUX\n" << FLUX << std::endl;
  //abort();
#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x+1; ++i) {
    mesh[i] = mesh[i] + dt * us[i];
  }
#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x; ++i) {
    Con[i] = (Con_n[i]*(mesh_n[i+1]-mesh_n[i]) - dt*(FLUX[i+1] - FLUX[i])) / (mesh[i+1]-mesh[i]);
    Pri[i] = Con2Pri(Con[i], Gamma[i]);
  }
  //stage 2
  Reconstruction(Con, mesh, ReconL_Con, ReconR_Con, ReconL_Pri, ReconR_Pri);
  cal_flux_LF(ReconL_Con, ReconR_Con, ReconL_Pri, ReconR_Pri, FLUX, alpha);
  cal_us_roeav(ReconL_Pri, ReconR_Pri, us);
  std::cout << "stage2" << std::endl;
  std::cout << "Con\n" << Con << std::endl;
  std::cout << "Pri\n" << Pri << std::endl;
  std::cout << "ReconL_Con\n" << ReconL_Con << std::endl;
  std::cout << "ReconR_Con\n" << ReconR_Con << std::endl;
  std::cout << "ReconL_Pri\n" << ReconL_Pri << std::endl;
  std::cout << "ReconR_Pri\n" << ReconR_Pri << std::endl;
  std::cout << "FLUX\n" << FLUX << std::endl;

#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x+1; ++i) {
    mesh1[i] = 0.75*mesh_n[i] + 0.25*(mesh[i] + dt * us[i]);
  }
#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x; ++i) {
    Con[i] = ( 0.75*Con_n[i]*(mesh_n[i+1]-mesh_n[i])
        + 0.25*(Con[i]*(mesh[i+1]-mesh[i]) - dt*(FLUX[i+1] - FLUX[i])) )/ (mesh1[i+1]-mesh1[i]);
    Pri[i] = Con2Pri(Con[i], Gamma[i]);
  }
  //stage 3
  Reconstruction(Con, mesh, ReconL_Con, ReconR_Con, ReconL_Pri, ReconR_Pri);
  cal_flux_LF(ReconL_Con, ReconR_Con, ReconL_Pri, ReconR_Pri, FLUX, alpha);
  cal_us_roeav(ReconL_Pri, ReconR_Pri, us);
  std::cout << "stage3" << std::endl;
  std::cout << "Con\n" << Con << std::endl;
  std::cout << "Pri\n" << Pri << std::endl;
  std::cout << "ReconL_Con\n" << ReconL_Con << std::endl;
  std::cout << "ReconR_Con\n" << ReconR_Con << std::endl;
  std::cout << "ReconL_Pri\n" << ReconL_Pri << std::endl;
  std::cout << "ReconR_Pri\n" << ReconR_Pri << std::endl;
  std::cout << "FLUX\n" << FLUX << std::endl;

#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x+1; ++i) {
    mesh[i] = 1./3*mesh_n[i] + 2./3*(mesh1[i] + dt * us[i]);
  }
#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x; ++i) {
    Con[i] = ( 1./3*Con_n[i]*(mesh_n[i+1]-mesh_n[i])
        + 2./3*(Con[i]*(mesh1[i+1]-mesh1[i]) - dt*(FLUX[i+1] - FLUX[i])) )/ (mesh[i+1]-mesh[i]);
    Pri[i] = Con2Pri(Con[i], Gamma[i]);
  }

}

void Eulerian1D::SSP_RK_HLLC(Sol& Con, Sol& Pri, VEC& mesh, const double dt) {
  VEC mesh_n(mesh), mesh1(N_x+1);
  Sol Con_n(Con), Pri_n(Pri);
  //stage 1
  Reconstruction(Con_n, mesh, ReconL_Con, ReconR_Con, ReconL_Pri, ReconR_Pri);
  cal_flux_HLLC(ReconL_Con, ReconR_Con, ReconL_Pri, ReconR_Pri, FLUX, us);
#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x+1; ++i) {
    mesh[i] = mesh[i] + dt * us[i];
  }
#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x; ++i) {
    Con[i] = (Con_n[i]*(mesh_n[i+1]-mesh_n[i]) - dt*(FLUX[i+1] - FLUX[i])) / (mesh[i+1]-mesh[i]);
    Pri[i] = Con2Pri(Con[i], Gamma[i]);
  }
  //stage 2
  Reconstruction(Con, mesh, ReconL_Con, ReconR_Con, ReconL_Pri, ReconR_Pri);
  cal_flux_HLLC(ReconL_Con, ReconR_Con, ReconL_Pri, ReconR_Pri, FLUX, us);
#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x+1; ++i) {
    mesh1[i] = 0.75*mesh_n[i] + 0.25*(mesh[i] + dt * us[i]);
  }
#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x; ++i) {
    Con[i] = ( 0.75*Con_n[i]*(mesh_n[i+1]-mesh_n[i])
        + 0.25*(Con[i]*(mesh[i+1]-mesh[i]) - dt*(FLUX[i+1] - FLUX[i])) )/ (mesh1[i+1]-mesh1[i]);
    Pri[i] = Con2Pri(Con[i], Gamma[i]);
  }
  //stage 3
  Reconstruction(Con, mesh, ReconL_Con, ReconR_Con, ReconL_Pri, ReconR_Pri);
  cal_flux_HLLC(ReconL_Con, ReconR_Con, ReconL_Pri, ReconR_Pri, FLUX, us);
#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x+1; ++i) {
    mesh[i] = 1./3*mesh_n[i] + 2./3*(mesh1[i] + dt * us[i]);
  }
#pragma omp parallel for num_threads(Nthread)
  for(u_int i = 0; i < N_x; ++i) {
    Con[i] = ( 1./3*Con_n[i]*(mesh_n[i+1]-mesh_n[i])
        + 2./3*(Con[i]*(mesh1[i+1]-mesh1[i]) - dt*(FLUX[i+1] - FLUX[i])) )/ (mesh[i+1]-mesh[i]);
    Pri[i] = Con2Pri(Con[i], Gamma[i]);
  }
}

