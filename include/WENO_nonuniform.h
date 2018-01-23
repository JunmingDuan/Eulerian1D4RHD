#ifndef _WENO_NONUNIFORM_H
#define _WENO_NONUNIFORM_H

#include <cmath>

/**
 * @brief WENO_cell_L weno reconstruction by cell average on non-uniform mesh
 *
 * @param h1 cell {i-1} length
 * @param h2
 * @param h3 cell {i+1} length
 * @param u1 cell {i-1} average
 * @param u2
 * @param u3 cell {i+1} average
 *
 * @return point val at x_{i-1/2}
 */
double WENO_cell_L(double h1, double h2, double h3, double u1, double u2, double u3);

/**
 * @brief WENO_cell_R weno reconstruction by cell average on non-uniform mesh
 *
 * @param h1 cell {i-1} length
 * @param h2
 * @param h3 cell {i+1} length
 * @param u1 cell {i-1} average
 * @param u2
 * @param u3 cell {i+1} average
 *
 * @return point val at x_{i+1/2}
 */
double WENO_cell_R(double h1, double h2, double h3, double u1, double u2, double u3);

#endif

