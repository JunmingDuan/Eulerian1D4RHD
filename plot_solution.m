DAT = load('sol.dat');
x = DAT(:,1);
rho = DAT(:,2);
u = DAT(:,3);
p = DAT(:,4);
e = DAT(:,5);

%figure(1)
%plot(x, rho/10, '-o');
%figure(2)
%plot(x, u, '-o');
%figure(3)
%plot(x, p*3/40, '-o');
%figure(4)
%plot(x, e, '-o');

figure(1)
%plot(x, rho/10, '-o');
plot(x, rho, '-o');
axis([0.49,0.54,-10,130]);
figure(2)
plot(x, u, '-o');
axis([0.49,0.54,-1,1.2]);
figure(3)
%plot(x, p*3/40, '-o');
%plot(x, p/1e3, '-o');
plot(x, p, '-o');
axis([0.49,0.54,-25,400]);
figure(4)
%plot(x, e/1.5e3, '-o');
plot(x, e, '-o');
axis([0.49,0.54,-50,800]);


