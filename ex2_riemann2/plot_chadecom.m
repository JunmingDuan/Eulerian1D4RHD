DAT1 = load('sol.dat');
x1 = DAT1(:,1);
rho1 = DAT1(:,2);
u1 = DAT1(:,3);
p1 = DAT1(:,4);
e1 = DAT1(:,5);

EX1 = load('../exact_solution/2011GRPex4.2.dat');
x0 = EX1(:,1);
p0 = EX1(:,2);
rho0 = EX1(:,3);
u0 = EX1(:,4);
e0 = EX1(:,5);

figure(1)
plot(x1, rho1, 'or', x0, rho0, '-k');
figure(2)
plot(x1, u1, 'or', x0, u0, '-k');
figure(3)
plot(x1, p1, 'or', x0, p0, '-k');
figure(4)
plot(x1, e1, 'or', x0, e0, '-k');

