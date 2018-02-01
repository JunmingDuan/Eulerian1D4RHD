clear all;
DAT1 = load('ex5_HLLC_n10000_RK3_Cha_Eul.dat');
%DAT1 = load('sol.dat');
x1 = DAT1(:,1);
rho1 = DAT1(:,2);
u1 = DAT1(:,3);
p1 = DAT1(:,4);
e1 = DAT1(:,5);

DAT2 = load('ex5_HLLC_n4000_RK3_Cha_Eul.dat');
%DAT2 = load('sol.dat');
x2 = DAT2(:,1);
rho2 = DAT2(:,2);
u2 = DAT2(:,3);
p2 = DAT2(:,4);
e2 = DAT2(:,5);


figure(1)
plot(x1, rho1, 'or', x2, rho2, '-b', x1, rho1, '-k');
legend('Lag', 'Eul', 'Location', 'NorthEast');
axis([0.49, 0.54, -10, 125]);
figure(2)
plot(x1, u1, 'or', x2, u2, '-b', x1, u1, '-k');
legend('Lag', 'Eul', 'Location', 'NorthWest');
axis([0.49, 0.54, -0.9, 1.1]);
figure(3)
plot(x1, p1, 'or', x2, p2, '-b', x1, p1, '-k');
legend('Lag', 'Eul', 'Location', 'NorthEast');
axis([0.49, 0.54, -25, 400]);
figure(4)
plot(x1, e1, 'or', x2, e2, '-b', x1, e1, '-k');
legend('Lag', 'Eul', 'Location', 'NorthEast');
axis([0.49, 0.54, -50, 800]);

