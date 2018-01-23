#ifndef _PARA_H
#define _PARA_H

extern int Nthread;
//0 for piecewise constant(not reconstruct), 1 for conservative variables, 2 for primitive
//variables, 3 for characteristic variables
extern int is_RECON;
//1 for outflow, -1 for reflextion
extern int BD;

#endif //_PARA_H

