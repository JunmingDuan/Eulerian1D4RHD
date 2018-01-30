#ifndef _ENO_H_
#define _ENO_H_
/**
 * @file ENO.h
 * @brief ENO reconstruction for 1D lagrangian scheme
 * @author Duan Junming, duanjm@pku.edu.cn
 * @version 1.0
 * @date 2017-07-16
 */
#include <vector>
#include <cmath>
#include <cstdlib>

/**
 * @brief ENO2
 *
 * @param h1 I_{i-1} length
 * @param h2
 * @param h3 I_{i+1} length
 * @param u1 I_{i-1} cell average
 * @param u2
 * @param u3 I_{i+1} cell average
 * @param point Positivity preserving on these points, these points is scale in [0, 1].
 *
 * @return u_{i-1/2}^+ and u_{i+1/2}^-
 */
template <class T>
void ENO2(double h1, double h2, double h3, const T& u1, const T& u2, const T& u3, const std::vector<double>& point,
    T& ul, T& ur);

double ENO2_CELL_L(double h1, double h2, double h3, const double u1, const double u2, const double u3);
double ENO2_CELL_R(double h1, double h2, double h3, const double u1, const double u2, const double u3);

void ENO2(double h1, double h2, double h3, const double u1, const double u2, const double u3,
    double& ul, double& ur);

#endif

