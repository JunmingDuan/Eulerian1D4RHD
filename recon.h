#ifndef __RECON_H
#define __RECON_H

#include "vvector.h"
#include "weno.cpp"

template<class T>
class recon {
  private:
		typedef vvector<vvector<T> > Sol;
  public:
    recon() {}

    //bm = -1, reflex; bm = 1. period;
    static void recon_x(const Sol& u, Sol& recon_ul, Sol& recon_ur, const int bm) {
      int nx = u.size();
      int ny = u[0].size();
      recon_ul.resize(nx);
      recon_ur.resize(nx);
      for(int i = 0; i < nx; ++i) {
        recon_ul[i].resize(ny);
        recon_ur[i].resize(ny);
      }
      for(int j = 0; j < ny; ++j) {
        for(int i = 0; i < nx; ++i) {
          if(i < nx-1 && i > 0) recon_ul[i][j] = WENO3CellR(u[i-1][j], u[i][j], u[i+1][j]);
          else if(i == nx-1 && bm == 1) recon_ul[i][j] = WENO3CellR(u[i-1][j], u[i][j], u[0][j]);
          else if(i == 0 && bm == 1) recon_ul[i][j] = WENO3CellR(u[nx-1][j], u[i][j], u[i+1][j]);
          else std::cout << "Wrong!" << std::endl;
          if(i > 0 && i < nx-1) recon_ur[i][j] = WENO3CellL(u[i-1][j], u[i][j], u[i+1][j]);
          else if(i == nx-1 && bm == 1) recon_ur[i][j] = WENO3CellL(u[i-1][j], u[i][j], u[0][j]);
          else if(i == 0 && bm == 1) recon_ur[i][j] = WENO3CellL(u[nx-1][j], u[i][j], u[i+1][j]);
          else std::cout << "Wrong!" << std::endl;
        }
      }
    }

    static void recon_y(const Sol& u, Sol& recon_ul, Sol& recon_ur, const int bm) {
      int nx = u.size();
      int ny = u[0].size();
      recon_ul.resize(nx);
      recon_ur.resize(nx);
      for(int i = 0; i < nx; ++i) {
        recon_ul[i].resize(ny);
        recon_ur[i].resize(ny);
      }
      for(int i = 0; i < nx; ++i) {
        for(int j = 0; j < ny; ++j) {
          if(j < ny-1 && j > 0) recon_ul[i][j] = WENO3CellR(u[i][j-1], u[i][j], u[i][j+1]);
          else if(j == ny-1 && bm == 1) recon_ul[i][j] = WENO3CellR(u[i][j-1], u[i][j], u[i][0]);
          else if(j == 0 && bm == 1) recon_ul[i][j] = WENO3CellR(u[i][ny-1], u[i][j], u[i][j+1]);
          else std::cout << "Wrong!" << std::endl;
          if(j > 0 && j < ny-1) recon_ur[i][j] = WENO3CellL(u[i][j-1], u[i][j], u[i][j+1]);
          else if(j == ny-1 && bm == 1) recon_ur[i][j] = WENO3CellL(u[i][j-1], u[i][j], u[i][0]);
          else if(j == 0 && bm == 1) recon_ur[i][j] = WENO3CellL(u[i][ny-1], u[i][j], u[i][j+1]);
          else std::cout << "Wrong!" << std::endl;
        }
      }

    }

};

#endif

