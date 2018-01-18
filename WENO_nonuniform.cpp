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
double WENO_cell_L(double h1, double h2, double h3, double u1, double u2, double u3) {
  double eps = 1e-6;
  double h22 = pow(h2, 2);
  double t1 = pow(u1-u2, 2);
  double t2 = pow(u2-u3, 2);
  double t3 = pow(h1+h2, 2);
  double t4 = pow(h2+h3, 2);
  return
  (h1*(2*h2*u2 - h2*u3 + h3*u2))/((h2 + h3)*((h2 + h3)/(pow(eps + (4*h22*t1)/t3, 2)*(h1 + h2 + h3))
        + h1/(pow(eps + (4*h22*t2)/t4, 2)*(h1 + h2 + h3)))*pow(eps + (4*h22*t2)/t4, 2)*(h1 + h2 + h3))
  + ((h2 + h3)*(h1*u2 + h2*u1))/((h1 + h2)*((h2 + h3)/(pow(eps + (4*h22*t1)/t3, 2)*(h1 + h2 + h3))
        + h1/(pow(eps + (4*h22*t2)/t4, 2)*(h1 + h2 + h3)))*pow(eps + (4*h22*t1)/t3, 2)*(h1 + h2 + h3));
}

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
double WENO_cell_R(double h1, double h2, double h3, double u1, double u2, double u3) {
  double eps = 1e-6;
  double h22 = pow(h2, 2);
  double t1 = pow(u1-u2, 2);
  double t2 = pow(u2-u3, 2);
  double t3 = pow(h1+h2, 2);
  double t4 = pow(h2+h3, 2);
  return
(h3*(h1*u2 - h2*u1 + 2*h2*u2))/((h1 + h2)*((h1 + h2)/(pow(eps + (4*h22*t2)/t4, 2)*(h1 + h2 + h3))
      + h3/(pow(eps + (4*h22*t1)/t3, 2)*(h1 + h2 + h3)))*pow(eps + (4*h22*t1)/t3, 2)*(h1 + h2 + h3))
+ ((h1 + h2)*(h2*u3 + h3*u2))/((h2 + h3)*((h1 + h2)/(pow(eps + (4*h22*t2)/t4, 2)*(h1 + h2 + h3))
      + h3/(pow(eps + (4*h22*t1)/t3, 2)*(h1 + h2 + h3)))*pow(eps + (4*h22*t2)/t4, 2)*(h1 + h2 + h3));
}

