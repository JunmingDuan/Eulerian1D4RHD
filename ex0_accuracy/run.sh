#!/bin/sh

make -j;
for n in 25 50 100 200 400 800 1600;
do
  ./main $n 0.4 2;
  #mv sol.dat ex0_LF_n${n}_ENO2.dat;
  #cat err.dat >> LF_err_ENO2.dat;
  mv sol.dat ex0_HLLC_n${n}_WENO3.dat;
  cat err.dat >> HLLC_err_WENO3.dat;
done


