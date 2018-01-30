#include "Eulerian1D.h"

void Eulerian1D::ROE_AV_MAT(const bU& PRIL, const bU& PRIR, const double GAMMAL, const double GAMMAR,
    MAT& R, MAT& L) {
    double hl, hr, kl, kr, gl, gr, ul, ur;
    double v0, v1, v2;
    hl = 1 + GAMMAL/(GAMMAL-1)*PRIL[2]/PRIL[0];
    hr = 1 + GAMMAR/(GAMMAR-1)*PRIR[2]/PRIR[0];
    kl = sqrt(PRIL[0]*hl);
    kr = sqrt(PRIR[0]*hr);
    ul = PRIL[1];
    ur = PRIR[1];
    gl = 1./sqrt(1-ul*ul);
    gr = 1./sqrt(1-ur*ur);
    v0 = (kl*gl + kr*gr)/(kl+kr);
    v1 = (kl*gl*ul + kr*gr*ur)/(kl+kr);
    v2 = (kl*PRIL[2]/PRIL[0]/hl + kr*PRIR[2]/PRIR[0]/hr) / (kl+kr);
    double cm = 1 - GAMMAL*v2/(GAMMAL-1);
    double cp = 1 + GAMMAL*v2/(GAMMAL-1);
    double va = - v0*v0 + v1*v1;
    double s2 = 0.5*GAMMAL*v2*(1-va) - 0.5*(GAMMAL-1)*(1+va);
    double s = sqrt(s2);
    double e = -va;
    double y = sqrt((1-GAMMAL*v2)*e + s2);

    {//right characteristic matrix, R
      R[0][0] = cm;          R[0][1] = cm + s2/(GAMMAL-1); R[0][2] = cm;
      R[1][0] = v1 - s/y*v0; R[1][1] = v1;                 R[1][2] = v1 + s/y*v0;
      R[2][0] = v0 - s/y*v1; R[2][1] = v0;                 R[2][2] = v0 + s/y*v1;
    }
    {//left characteristic matrix, L
      L[0][0] = (GAMMAL-1)*e;    L[0][1] = - s2*v1+s*y*v0 + (GAMMAL-1)*e*cp*v1; L[0][2] = s2*v0 - s*y*v1 - (GAMMAL-1)*e*cp*v0;
      L[1][0] = -2*(GAMMAL-1)*e; L[1][1] = 4*s2*v1 - 2*(GAMMAL-1)*e*cp*v1;      L[1][2] = -4*s2*v0 + 2*(GAMMAL-1)*e*cp*v0;
      L[2][0] = (GAMMAL-1)*e;    L[2][1] = - s2*v1-s*y*v0 + (GAMMAL-1)*e*cp*v1; L[2][2] = s2*v0 + s*y*v1 - (GAMMAL-1)*e*cp*v0;
      L *= -0.5/e/s2;
    }
    //R[0][0] = 1; R[0][1] = 0; R[0][2] = 0;
    //R[1][0] = 0; R[1][1] = 1; R[1][2] = 0;
    //R[2][0] = 0; R[2][1] = 0; R[2][2] = 1;
    //L[0][0] = 1; L[0][1] = 0; L[0][2] = 0;
    //L[1][0] = 0; L[1][1] = 1; L[1][2] = 0;
    //L[2][0] = 0; L[2][1] = 0; L[2][2] = 1;
}

bU Eulerian1D::multiply(const bU& x, MAT& M) {
  bU y;
  y[0] = M[0][0]*x[0] + M[0][1]*x[1] + M[0][2]*x[2];
  y[1] = M[1][0]*x[0] + M[1][1]*x[1] + M[1][2]*x[2];
  y[2] = M[2][0]*x[0] + M[2][1]*x[1] + M[2][2]*x[2];
  return y;
}

/*void X_EigenMatrix(double V[], double EigR1vecx[][DOS], double EigL1vecx[][DOS])*/
//{
  //const double rho = V[0];
  //const double u1 = V[1];
  //const double u2 = V[2];
  //const double p = V[3];
  //const double rhoh = rho + GAMMA/(GAMMA-1.) * p;
  //const double h = rhoh / rho;
  //const double cs2 = GAMMA*p / rhoh;
  //const double cs = sqrt(cs2);
  //const double u12 = u1 * u1;
  //const double u22 = u2 * u2;
  //const double unorm2 = u12 + u22;
  //const double gamma2 = 1.0 / ( 1.0 - u12 - u22 );
  //const double gamma = sqrt(gamma2);
  //const double G2 = 1. / (GAMMA-1.);
  //double Lambda[DOS];
  //double wk;
  //int k;

  //Lambda[0] = u1*(1-cs2) - cs * sqrt( (1.-unorm2)*(1.-u12-u22*cs2) );
  //Lambda[0] = Lambda[0] / ( 1-unorm2*cs2 );
  //Lambda[1] = u1;
  //Lambda[2] = u1;
  //Lambda[3] = u1*(1-cs2) + cs * sqrt( (1.-unorm2)*(1.-u12-u22*cs2) );
  //Lambda[3] = Lambda[3] / ( 1-unorm2*cs2 );

  //EigR1vecx[0][0] = rho * gamma * (u1 * Lambda[0]-1.0);
  //EigR1vecx[0][1] = 1.0;
  //EigR1vecx[0][2] = u2 * gamma;
  //EigR1vecx[0][3] = rho * gamma * (u1 * Lambda[3]-1.0);

  //EigR1vecx[1][0] = gamma2 * rhoh * Lambda[0] * (u12 - 1.0);
  //EigR1vecx[1][1] = u1 * gamma;
  //EigR1vecx[1][2] = 2 * h * u1 * u2 * gamma2;
  //EigR1vecx[1][3] = gamma2 * rhoh * Lambda[3] * (u12 - 1.0);

  //EigR1vecx[2][0] = gamma2 * rhoh * u2 * ( u1*Lambda[0]-1.0 );
  //EigR1vecx[2][1] = u2 * gamma;
  //EigR1vecx[2][2] = h * ( 1.0 + 2*u22*gamma2 );
  //EigR1vecx[2][3] = gamma2 * rhoh * u2 * ( u1*Lambda[3]-1.0 );

  //EigR1vecx[3][0] = gamma2 * rhoh * ( u12 - 1.0 );
  //EigR1vecx[3][1] = gamma;
  //EigR1vecx[3][2] = 2 * h * u2 * gamma2;
  //EigR1vecx[3][3] = gamma2 * rhoh * ( u12 - 1.0 );


  ////Left Eigen Matrix
  //k = 0;
  //wk = 2 * rhoh * G2 * ( u1-Lambda[k] ) * ( u1-Lambda[k] ) * gamma2 * ( u12 + u22*cs2 -1.0 );
  //cs2wk = ( 1.0 + cs2 * G2 ) / wk;
  //EigL1vecx[k][0] = wk * ( u1*Lambda[k]-1.0 ) / ( gamma * ( 1.0 + cs2*G2 ));
  //EigL1vecx[k][1] = Lambda[k] * ( (unorm2+G2)/(1.0+cs2*G2) - u22 ) - u1 * ( 1.0+G2 ) / ( 1.0+cs2*G2 );
  //EigL1vecx[k][1] = wk * EigL1vecx[k][1];
  //EigL1vecxEigL1vecx[k][2] = wk * u2 * ( u1 * Lambda[k] - 1.0 );
  //EigL1vecx[k][3] = u22 + 1.0 / ( gamma2*(1+cs2*G2) ) - u1 * ( Lambda[k]-u1 ) * ( 1.0+G2 ) / (1.0+cs2*G2);
  //EigL1vecx[k][3] = wk * EigL1vecx[k][3];

  //EigL1vecxEigL1vecx[1][0] = 1.0 / (cs2*G2);
  //EigL1vecx[1][1] = u1 * ( 1.0+u22*gamma2 ) / ( (1.0-u12)*gamma*h*cs2*G2 );
  //EigL1vecx[1][2] = u2 * gamma / ( G2*h*cs2 );
  //EigL1vecx[1][3] = ( 1.0+u22*gamma2 ) / ( (u12-1.0)*gamma*G2*h*cs2 );

  //EigL1vecx[2][0] = 0.0;
  //EigL1vecx[2][1] = u1*u2 / ( h*(1.0-u12) );
  //EigL1vecx[2][2] = 1.0 / h;
  //EigL1vecx[2][3] = u2 / ( h*(u12-1.0) );

  //k = 3;
  //wk = 2 * rhoh / (GAMMA-1.) * ( u1-Lambda[k] ) * ( u1-Lambda[k] ) * gamma2 * ( u12 + u22*cs2 -1.0 );
  //wk = ( 1.0 + cs2 * G2 ) / wk;
  //EigL1vecx[k][0] = wk * ( u1*Lambda[k]-1.0 ) / (gamma * ( 1.0 + cs2*G2 ));
  //EigL1vecx[k][1] = Lambda[k] * ( (unorm2+G2)/(1.0+cs2*G2) - u22 ) - u1 * ( 1.0+G2 ) / ( 1.0+cs2*G2 );
  //EigL1vecx[k][1] = wk * EigL1vecx[k][1];
  //EigL1vecx[k][2] = wk * u2 * ( u1 * Lambda[k] - 1.0 );
  //EigL1vecx[k][3] = u22 + 1.0 / ( gamma2*(1+cs2*G2) ) - u1 * ( Lambda[k]-u1 ) * ( 1.0+G2 ) / (1.0+cs2*G2);
  //EigL1vecx[k][3] = wk * EigL1vecx[k][3];
//}

//void Y_EigenMatrix(double V[], double REM[][DOS], double LEM[][DOS])
//{
  //const double rho = V[0];
  //const double u1 = V[2];
  //const double u2 = -V[1];
  //const double p = V[3];
  //const double rhoh = rho + GAMMA/(GAMMA-1.) * p;
  //const double h = rhoh / rho;
  //const double cs2 = GAMMA*p / rhoh;
  //const double cs = sqrt(cs2);
  //const double u12 = u1 * u1;
  //const double u22 = u2 * u2;
  //const double unorm2 = u12 + u22;
  //const double gamma2 = 1.0 / ( 1.0 - u12 - u22 );
  //const double gamma = sqrt(gamma2);
  //const double G2 = 1. / (GAMMA-1.);
  //double nx = 0.;
  //double ny = 1.;
  //int k;
  //Eig - information
  //informationdouble Lambda[DOS];
  //DOSdouble wk;
  //wkdouble EigR1vecx[DOS][DOS];
  //DOSdouble EigL1vecx[DOS][DOS];
  
  //DOSLambda[0] = u1*(1-cs2) - cs * sqrt( (1.-unorm2)*(1.-u12-u22*cs2) );
  //cs2Lambda[0] = Lambda[0] / ( 1-unorm2*cs2 );
  //cs2Lambda[1] = u1;
  //u1Lambda[2] = u1;
  //u1Lambda[3] = u1*(1-cs2) + cs * sqrt( (1.-unorm2)*(1.-u12-u22*cs2) );
  //cs2Lambda[3] = Lambda[3] / ( 1-unorm2*cs2 );
  
  //cs2EigR1vecx[0][0] = rho * gamma * (u1 * Lambda[0]-1.0);
  //LambdaEigR1vecx[0][1] = 1.0;
  //LambdaEigR1vecxEigR1vecx[0][2] = u2 * gamma;
  //gammaEigR1vecx[0][3] = rho * gamma * (u1 * Lambda[3]-1.0);
  
  //LambdaEigR1vecx[1][0] = gamma2 * rhoh * Lambda[0] * (u12 - 1.0);
  //u12EigR1vecx[1][1] = u1 * gamma;
  //gammaEigR1vecx[1][2] = 2 * h * u1 * u2 * gamma2;
  //gamma2EigR1vecx[1][3] = gamma2 * rhoh * Lambda[3] * (u12 - 1.0);
  
  //u12EigR1vecx[2][0] = gamma2 * rhoh * u2 * ( u1*Lambda[0]-1.0 );
  //LambdaEigR1vecx[2][1] = u2 * gamma;
  //gammaEigR1vecx[2][2] = h * ( 1.0 + 2*u22*gamma2 );
  //gamma2EigR1vecx[2][3] = gamma2 * rhoh * u2 * ( u1*Lambda[3]-1.0 );
  
  //LambdaEigR1vecx[3][0] = gamma2 * rhoh * ( u12 - 1.0 );
  //u12EigR1vecx[3][1] = gamma;
  //gammaEigR1vecx[3][2] = 2 * h * u2 * gamma2;
  //gamma2EigR1vecx[3][3] = gamma2 * rhoh * ( u12 - 1.0 );
  
  //u12k = 0;
  //u12kwk = 2 * rhoh * G2 * ( u1-Lambda[k] ) * ( u1-Lambda[k] ) * gamma2 * (
  //u12 + u22*cs2 -1.0 );
  //cs2wk = ( 1.0 + cs2 * G2 ) / wk;
  //wkEigL1vecx[k][0] = wk * ( u1*Lambda[k]-1.0 ) / ( gamma * ( 1.0 + cs2*G2 )
  //);
  //G2EigL1vecx[k][1] = Lambda[k] * ( (unorm2+G2)/(1.0+cs2*G2) - u22 ) - u1 * (
  //1.0+G2 ) / ( 1.0+cs2*G2 );
  //G2EigL1vecx[k][1] = wk * EigL1vecx[k][1];
  //EigL1vecxEigL1vecx[k][2] = wk * u2 * ( u1 * Lambda[k] - 1.0 );
  //LambdaEigL1vecx[k][3] = u22 + 1.0 / ( gamma2*(1+cs2*G2) ) - u1 * (
  //Lambda[k]-u1 ) * ( 1.0+G2 ) / (1.0+cs2*G2);
  //G2EigL1vecx[k][3] = wk * EigL1vecx[k][3];
  
  //EigL1vecxEigL1vecx[1][0] = 1.0 / (cs2*G2);
  //G2EigL1vecx[1][1] = u1 * ( 1.0+u22*gamma2 ) / ( (1.0-u12)*gamma*h*cs2*G2 );
  //G2EigL1vecx[1][2] = u2 * gamma / ( G2*h*cs2 );
  //cs2EigL1vecx[1][3] = ( 1.0+u22*gamma2 ) / ( (u12-1.0)*gamma*G2*h*cs2 );
  
  //cs2EigL1vecx[2][0] = 0.0;
  //cs2EigL1vecxEigL1vecx[2][1] = u1*u2 / ( h*(1.0-u12) );
  //u12EigL1vecx[2][2] = 1.0 / h;
  //u12EigL1vecxEigL1vecx[2][3] = u2 / ( h*(u12-1.0) );
  
  //u12k = 3;
  //u12kwk = 2 * rhoh * G2 * ( u1-Lambda[k] ) * ( u1-Lambda[k] ) * gamma2 * (
  //u12 + u22*cs2 -1.0 );
  //cs2wk = ( 1.0 + cs2 * G2 ) / wk;
  //wkEigL1vecx[k][0] = wk * ( u1*Lambda[k]-1.0 ) / (gamma * ( 1.0 + cs2*G2 ));
  //G2EigL1vecx[k][1] = Lambda[k] * ( (unorm2+G2)/(1.0+cs2*G2) - u22 ) - u1 * (
  //1.0+G2 ) / ( 1.0+cs2*G2 );
  //G2EigL1vecx[k][1] = wk * EigL1vecx[k][1];
  //EigL1vecxEigL1vecx[k][2] = wk * u2 * ( u1 * Lambda[k] - 1.0 );
  //LambdaEigL1vecx[k][3] = u22 + 1.0 / ( gamma2*(1+cs2*G2) ) - u1 * (
  //Lambda[k]-u1 ) * ( 1.0+G2 ) / (1.0+cs2*G2);
  //G2EigL1vecx[k][3] = wk * EigL1vecx[k][3];
  
  //EigL1vecx[>计算沿第ie条边外法方向的特征矩阵，在这里利用旋转不变性
  //(n_x,n_y)=(cos(theta),sin(theta))
  //thetaR(theta) = T(theta)^-1*R
  //thetaL(theta) = L * T(theta)
  //thetaT(theta) = ( 1     0            0      0 )     T(theta)^-1 = ( 1     0
  //0       0 )
  //theta( 0 cos(theta)   sin(theta) 0 )                   ( 0 cos(theta)
  //-sin(theta) 0 )
  //theta( 0 -sin(theta)) cos(theta) 0 )                   ( 0 sin(theta))
  //cos(theta) 0 )
  //theta( 0     0            0      1 )                   ( 0     0
  //0      1 )
  //theta*/
  //theta//EigLvecxRot
  //EigLvecxRotfor(int col = 0; col < DOS; col++){
  //DOSswitch(col){
  //DOScase 0: case 3:
  //colfor(int row = 0; row < DOS; row++){
  //colforLEM[row][col] = EigL1vecx[row][col];
  //EigL1vecx}
  //rowbreak;
  //colcase 1:
  //EigL1vecxfor(int k = 0; k < DOS; k++){
  //colcaseLEM[k][col] = EigL1vecx[k][1]*nx - EigL1vecx[k][2]*ny;
  //nx}
  //EigL1vecxbreak;
  //nycase 2:
  //nxfor(int k = 0; k < DOS; k++){
  //nycaseLEM[k][col] = EigL1vecx[k][1]*ny + EigL1vecx[k][2]*nx;
  //ny}
  //EigL1vecxbreak;
  //nxdefault:
  //nycout<<"The index of collum is out of the randge!"<<endl;
  //theexit(1);
  //randgebreak;
  //theexit}
  //theexit}
  
  //theexit//EigRvecxRot
  //EigRvecxRotfor(int row = 0; row < DOS; row++){
  //DOSswitch(row){
  //DOScase 0: case 3:
  //rowfor(int col = 0; col < DOS; col++){
  //rowforREM[row][col] = EigR1vecx[row][col];
  //EigR1vecx}
  //rowbreak;
  //colcase 1:
  //EigR1vecxfor(int col = 0; col < DOS; col++){
  //EigR1vecxforREM[row][col] = EigR1vecx[1][col]*nx - EigR1vecx[2][col]*ny;
  //EigR1vecx}
  //colbreak;
  //nycase 2:
  //EigR1vecxfor(int col = 0; col < DOS; col++){
  //EigR1vecxforREM[row][col] = EigR1vecx[1][col]*ny + EigR1vecx[2][col]*nx;
  //EigR1vecx}
  //colbreak;
  //nxdefault:
  //EigR1vecxcout<<"The index of collum is out of the randge!"<<endl;
  //theexit(1);
  //randgebreak;
  //theexit}
  //theexit}
  
  //}
  //}}}}}}}})]))]])])]]]]]]]}
