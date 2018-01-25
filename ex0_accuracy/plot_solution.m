clear all;
t_end = 2;
x = linspace(0, 1, 1e3);
y = 1+0.2*sin(2*pi*(x-0.99*t_end));

%DAT = load('ex0_LF_n800_WENO3.dat');
DAT = load('ex0_HLLC_n800_WENO3.dat');
x1 = DAT(:,1);
rho1 = DAT(:,2);

figure(1)
plot(x1, rho1, 'or', x, y, '-k');
legend('Numer', 'exact');

%err = load('LF_err_WENO3.dat');
err = load('HLLC_err_WENO3.dat');
disp('    L1        L2        Linf');
err_order = [
-log(err(2:end,2)./err(1:end-1,2))./log(err(2:end,1)./err(1:end-1,1)), ...
-log(err(2:end,3)./err(1:end-1,3))./log(err(2:end,1)./err(1:end-1,1)), ...
-log(err(2:end,4)./err(1:end-1,4))./log(err(2:end,1)./err(1:end-1,1))
];
disp(err_order);

