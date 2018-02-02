%two blast wave interaction
%final time = 0.43
%rho = 1, u = 0
%x < 0.1: p = 1e3
%x < 0.9: p = 1e-2
%x < 1.0: p = 1e2
t = 0.43;

R1_rho = 1.000e+0;
R1_u   = 0.000e-0;
R1_p   = 1.000e+3;
R1_c   = 6.323e-1;

R3_rho = 4.910e-2;
R3_u   = 9.570e-1;
R3_p   = 1.471e+1;
R3_c   = 6.321e-1;

R4_rho = 1.439e+1;
R4_u   = 9.570e-1;
R4_p   = 1.471e+1;
R4_c   = 5.591e-1;

R6_rho = 9.720e+0;
R6_u   = -8.820e-1;
R6_p   = 4.639e+0;
R6_c   = 5.002e-1;

R7_rho = 1.120e-1;
R7_u   = -8.820e-1;
R7_p   = 4.639e+0;
R7_c   = 6.303e-1;

R9_rho = 1.000e+0;
R9_u   = 0.000e-0;
R9_p   = 1.000e+2;
R9_c   = 6.316e-1;

C1_rho = 1.044e+2;
C1_u   = 4.560e-1;
C1_p   = 3.698e+2;
C1_c   = 6.084e-1;

C2_rho = 1.173e+2;
C2_u   = 4.560e-1;
C2_p   = 3.698e+2;
C2_c   = 6.056e-1;

x1 = 0.1 + 0.6324*t;
x2 = 0.1 + 0.8222*t;
x3 = 0.1 + 0.9570*t;
x4 = 0.1 + 0.9776*t;
x5 = 0.9 - 0.9274*t;
x6 = 0.9 - 0.8820*t;
x7 = 0.9 - 0.5668*t;
x8 = 0.9 - 0.6315*t;
x41 = 0.5106 + 0.088*(t - 0.42);
x42 = 0.5106 + 0.456*(t - 0.42);
x43 = 0.5106 + 0.703*(t - 0.42);
x1
x2
x3
x41
x42
x43
x6
x7
x8

x = [ linspace(0, x1, 100), linspace(x1, x2, 100), linspace(x2, x3, 100), linspace(x3, x4, 100), ...
linspace(x41,x42,100), linspace(x42,x43,100), linspace(x43,x6, 100), linspace(x6, x7, 100), ...
linspace(x7, x8, 100), linspace(x8, 1,  100)];
y = [ ones(1,100)*R1_rho, ones(1,100)*0, ones(1,100)*R3_rho, ones(1,100)*R4_rho, ...
ones(1,100)*C1_rho, ones(1,100)*C2_rho, ones(1,100)*R6_rho, ones(1,100)*R7_rho, ...
ones(1,100)*0, ones(1,100)*R9_rho];
plot(x, y);

