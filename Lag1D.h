/**
 * @file SCL1D.h
 * @brief 1D Lagrangian scheme for relativistic Euler equations
 * @author Duan Junming, duanjm@pku.edu.cn
 * @version 1.0
 * @date 2018-01-15
 */
#ifndef SCL1D_H
#define SCL1D_H

#include <iostream>
#include <assert.h>
#include <cstdlib>
#include "vvector.h"

typedef mvector<double,3> bU; //(rho, u, E), cell average

double EOS(const bU& u, const double gamma) {
  double p;
  p = u[0]*(u[2] - 0.5*u[1]*u[1])*(gamma-1.0);
  return p;
}

class GHOST {
  public:
    bU u;
    double gamma;
};

template<class T>
class SCL1D {
	private:
		typedef T (*Lfun)(const double t, const double x, double& gamma);
		typedef T (*NLfun)(T u, double t, double x);
		typedef vvector<T> Sol;
		u_int N_x;
		double t_start;
		double t_end;
		double x_start;
		double x_end;
		Lfun initial;
		double CFL;
		Sol NumSol;//数值解
    vvector<double> mi;
    vvector<double> mesh;
    vvector<double> gamma;
    vvector<double> us;
    vvector<double> A;
    Sol flux;
    Sol source;
    GHOST ghostl, ghostr;

	public:
		SCL1D(u_int Nx, double t_start, double t_end, double x_start, double x_end,
        Lfun initial, double CFL = 0.5)
      : N_x(Nx), t_start(t_start), t_end(t_end), x_start(x_start), x_end(x_end),
      initial(initial), CFL(CFL) {
        NumSol.resize(N_x);
        mi.resize(N_x);
        mesh.assign(N_x+1, 0);
        gamma.assign(N_x, 0);
        us.assign(N_x+1, 0);
        A.resize(N_x);
        flux.resize(N_x);
        source.resize(N_x);
        double h = (x_end - x_start) / N_x;
        for(u_int i = 0; i < N_x+1; ++i) {
          mesh[i] = h*i;
        }
      }

  private:
    void Initial() {
      for(u_int j = 0; j != N_x; ++j) {
        NumSol[j] = initial(t_start, 0.5*(mesh[j]+mesh[j+1]), gamma[j]);
        mi[j] = NumSol[j][0] * 0.5*(pow(mesh[j+1], 2) - pow(mesh[j], 2));
        A[j] = mesh[j+1] - mesh[j];
      }
    }

    void InfiniteBD() {
      ghostl.u = NumSol[0];
      ghostl.u[1] = -ghostl.u[1];
      double p = 1e-13;
      ghostr.u[0] = 1;
      ghostr.u[1] = -1;
      ghostr.u[2] = p/(5./3-1.0)+0.5;
      ghostl.gamma = gamma[0];
      ghostr.gamma = gamma[N_x-1];
    }

    double t_step(const double CFL) {
      double eps(1e-2);
      double dt1(eps);
      for(u_int i = 0; i < N_x; ++i) {
        double ai = sqrt(gamma[i]*EOS(NumSol[i], gamma[i])/NumSol[i][0]);
        double tmp = CFL*(mesh[i+1]-mesh[i])/ai;
        if(tmp < dt1) dt1 = tmp;
      }
      double dt2(eps);
      for(u_int i = 0; i < N_x; ++i) {
        double tmp = 0.4*(mesh[i+1]-mesh[i])/std::max(fabs(us[i]), fabs(us[i+1]));
        if(tmp < dt2) dt2 = tmp;
      }
      double sigma(1);
      for(u_int i = 0; i < N_x; ++i) {
        double ei = NumSol[i][2] - 0.5*pow(NumSol[i][1],2);
        double tmp = std::min(1.0, 0.5*NumSol[i][0]*ei/EOS(NumSol[i], gamma[i]));
        if(tmp < sigma) sigma = tmp;
      }
      double dt3(eps);
      for(u_int i = 0; i < N_x; ++i) {
        double tmp = sigma*mi[i]/NumSol[i][0]/fabs(mesh[i+1]*us[i+1]-mesh[i]*us[i]);
        if(tmp < dt3) dt3 = tmp;
      }
      double dt4(eps);
      for(u_int i = 0; i < N_x; ++i) {
        double z = NumSol[i][0]*sqrt(gamma[i]*EOS(NumSol[i], gamma[i])/NumSol[i][0]);
        double tmp = mi[i]/(mesh[i+1]*z+mesh[i]*z);
        if(tmp < dt4) dt4 = tmp;
      }
      //double dt5(1);
      //for(u_int i = 0; i < N_x; ++i) {
        //double pi = EOS(NumSol[i], gamma[i]);
        //double mu = (mesh[i+1] - mesh[i])*(pi -
        //double z = Numsol[i][0]*sqrt(gamma[i]/NumSol[i][0]);
        //double tmp = mi[i]/(mesh[i+1]*z+mesh[i]*z);
        //if(tmp < dt4) dt4 = tmp;
      //}
      std::cout << "dt1: " << dt1 << "; dt2: " << dt2 << "; sigma: " << sigma << "; dt3: " << dt3 << "; dt4: " << dt4 << std::endl;

      //return std::min(std::min(dt1, dt2), std::min(dt3, dt4));
      return std::min(dt1, dt2);
      //return 1e-2*CFL;
    }

    void rp_solver(const bU& UL, const bU& UR, double& us, double& ps, const double gammal, const double gammar) {
      double pl = EOS(UL, gammal);
      double pr = EOS(UR, gammar);
      double ul = UL[1];
      double ur = UR[1];
      double al = sqrt(gammal*pl/UL[0]);
      double ar = sqrt(gammar*pr/UR[0]);
      if(al != al) std::cout << "al wrong: " << pl << " " << UL[0] << std::endl;
      if(ar != ar) std::cout << "ar wrong: " << pr << " " << UR[0] << std::endl;
      double zl = UL[0] * al;
      double zr = UR[0] * ar;
      double sum = zl + zr;
      us = (zl*ul + zr*ur)/sum - (pr-pl)/sum;
      ps = (zl*pr + zr*pl)/sum - zl*zr*(ur-ul)/sum;
    }

    void cal_flux() {
      vvector<double> ps(N_x+1);
      for(u_int i = 1; i < N_x; ++i) {
        rp_solver(NumSol[i-1], NumSol[i], us[i], ps[i], gamma[i-1], gamma[i]);
      }
      //r=0不动
      rp_solver(ghostl.u, NumSol[0], us[0], ps[0], ghostl.gamma, gamma[0]);
      rp_solver(NumSol[N_x-1], ghostr.u, us[N_x], ps[N_x], gamma[N_x-1], ghostr.gamma);
      for(u_int i = 0; i < N_x; ++i) {
        //if(i != N_x-1) {
          flux[i][0] = 0;
          flux[i][1] = -(mesh[i+1]*ps[i+1] - mesh[i]*ps[i])/mi[i];
          flux[i][2] = -(mesh[i+1]*ps[i+1]*us[i+1] - mesh[i]*ps[i]*us[i])/mi[i];
        //}
        //else {
          //flux[i][0] = 0;
          //flux[i][1] = -(- mesh[i]*ps[i])/mi[i] + (mesh[i+1]-mesh[i])*EOS(NumSol[i],gamma[i])/mi[i];
          //flux[i][2] = -(- mesh[i]*ps[i]*us[i])/mi[i];
        //}
      }
    }

    void cal_source() {
      for(u_int i = 0; i < N_x; ++i) {
        source[i][0] = 0;
        source[i][1] = A[i]*EOS(NumSol[i],gamma[i])/mi[i];
        source[i][2] = 0;
      }
    }

    void move_mesh(double dt) {
      mesh += dt * us;
      for(u_int i = 0; i < N_x; ++i) {
        A[i] = (mesh[i+1]-mesh[i]);
        if(A[i] < 0) {
          std::cout << "mesh wrong! " << mesh[i] << " " << mesh[i+1] << std::endl;
          abort();
        }
      }
    }

    void modify_rho() {
      for(u_int i = 0; i < N_x; ++i) {
        NumSol[i][0] = mi[i]/( 0.5*(pow(mesh[i+1], 2) - pow(mesh[i], 2)) );
      }
    }

    void forward(double dt) {
      InfiniteBD();
      cal_flux();
      cal_source();
      move_mesh(dt);
      NumSol += dt * flux;
      NumSol += dt * source;
      modify_rho();
    }

  public:
		void Solve() {
      double t_now(t_start), dt;
			Initial();
			while(t_now < t_end - 1e-10) {
        if(t_now < 2.6e-2) dt = t_step(1e-6);
				else dt = t_step(CFL);
        if(dt + t_now > t_end) dt = t_end - t_now;

        forward(dt);

				t_now += dt;
        std::cout << "t: " << t_now << " , dt: " << dt << std::endl;
			}
      //std::cout << mi << std::endl;
      //std::cout << mesh << std::endl;
		}

		template<class Type> friend std::ostream& operator<<(std::ostream&,SCL1D<Type>&);
};

template<class T>
std::ostream& operator<<(std::ostream& os, SCL1D<T>& H){
	os.setf(std::ios::scientific);
	for(u_int i=0; i < H.N_x; i++) {
		os << 0.5*(H.mesh[i]+H.mesh[i+1]) << " " << H.NumSol[i] << "\n";
	}
	return os;
}

#endif

