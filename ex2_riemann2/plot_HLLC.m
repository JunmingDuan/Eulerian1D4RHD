%DAT1 = load('ex2_LF_n400_RK2_Con_Eul.dat');
DAT1 = load('sol.dat');
x1 = DAT1(:,1);
rho1 = DAT1(:,2);
u1 = DAT1(:,3);
p1 = DAT1(:,4);
e1 = DAT1(:,5);

DAT2 = load('ex2_HLLC_n400_RK3_Cha_Eul.dat');
%DAT2 = load('sol.dat');
x2 = DAT2(:,1);
rho2 = DAT2(:,2);
u2 = DAT2(:,3);
p2 = DAT2(:,4);
e2 = DAT2(:,5);

EX1 = load('../exact_solution/2011GRPex4.3.dat');
x0 = EX1(:,1);
p0 = EX1(:,2);
rho0 = EX1(:,3);
u0 = EX1(:,4);
e0 = EX1(:,5);

figure(1)
plot(x1, rho1/10, 'or', x2, rho2/10, '-b', x0, rho0/10, '-k');
legend('Recon', 'unRecon', 'exact');
figure(2)
plot(x1, u1, 'or', x2, u2, '-b', x0, u0, '-k');
legend('Recon', 'unRecon', 'exact');
figure(3)
plot(x1, p1/1e3, 'or', x2, p2/1e3, '-b', x0, p0/1e3, '-k');
legend('Recon', 'unRecon', 'exact');
figure(4)
plot(x1, e1/1.5e3, 'or', x2, e2/1.5e3, '-b', x0, e0/1.5e3, '-k');
legend('Recon', 'unRecon', 'exact');

