close all 
clear 
clc


prog = spotsosprog();

[prog,x] = prog.newIndeterminate('x',2);


q0 = [pi;0];
Q = diag([1 1]);
R = 1;
system = Pendulum(q0,Q,R);

f = system.Poly_f(x,-system.K*x);

V = .5*x'*system.S*x;
Vdot = diff(V,x)*f;

% found in VDP roa fixed v (SOSTest)
[prog, rho] = prog.newFree(1);
[prog, sig] = prog.newSOSPoly(monomials(x,2));

Vdot_sos = (x'*x)^3*(V - rho) - (1+sig)*Vdot;
prog = prog.withSOS(Vdot_sos);

spot_options = spotprog.defaultOptions;
spot_options.verbose = true;
solver = @spot_mosek;
sol = prog.minimize(-rho,solver,spot_options);



N = 30;
[X1,X2] = meshgrid(linspace(-5,5,N), linspace(-5,5,N));
X1DOT = X2;
X2DOT = -system.constants.b*X2-9.81*sin(X1)-(system.K(1).*X1+system.K(2).*X2); 
figure('Position',[0 0 1024 1024])
quiver(X1(:),X2(:),X1DOT(:),X2DOT(:),3.5)
xlim([-5 5])
ylim([-5 5])

hold on
VPLOT = reshape(dmsubs(V,x,[X1(:) X2(:)]'),size(X1));
[~,h] = contour(X1,X2,VPLOT,double(sol.eval(rho))*[1 1]);

set(h,'Color','Red','LineWidth',3)