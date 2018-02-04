syms r u p G rt ut pt rx ux px cs2
g = 1/sqrt(1-u^2);
h = 1 + p/r*G/(G-1);
D = r*g;
m = D*h*g*u;
E = D*h*g - p;

gt = g^3*u*ut;
ht = G/(G-1)*(pt*r - rt*p)/r^2;
Dt = rt*g + r*gt;
mt = Dt*h*g*u + D*ht*h*u + D*h*gt*u + D*h*g*ut;
Et = Dt*h*g + D*ht*g + D*h*gt - pt;

gx = g^3*u*ux;
hx = G/(G-1)*(px*r - rx*p)/r^2;
Dx = rx*g + r*gx;
mx = Dx*h*g*u + D*hx*h*u + D*h*gx*u + D*h*g*ux;
Ex = Dx*h*g + D*hx*g + D*h*gx - p;

eq1 = Dt + Dx*u + D*ux;
eq2 = mt + mx*u + m*ux + px;
eq3 = Et + mx;
eqns = [eq1==0, eq2==0, eq3==0];

sol = solve(eqns, [rt ut pt]);
%subs(sol.rt, G*p/r/(1+ p/r*G/(G-1)), cs2)
%subs(sol.ut, G*p/r/(1+ p/r*G/(G-1)), cs2)
%subs(sol.pt, G*p/r/(1+ p/r*G/(G-1)), cs2)
sol.rt = simplify(sol.rt);
sol.ut = simplify(sol.ut);
sol.pt = simplify(sol.pt);
pretty(diff(sol.ut, rx))
%pretty(subs(diff(sol.rt, rx), 1/(1-u^2)^(1/2), g));
%pretty(subs(sol.ut, 1/(1-u^2)^(1/2), g))
%pretty(subs(sol.pt, 1/(1-u^2)^(1/2), g))

