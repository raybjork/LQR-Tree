close all 
clear 
clc


prog = spotsosprog();

[prog,x] = prog.newIndeterminate('x',2);
%[prog,s] = prog.newIndeterminate('s',1);
%[prog,c] = prog.newIndeterminate('c',1);

%xhat = [x;s;c];
q0 = [pi;0];
Q = diag([1 1]);
R = 1;
system = Pendulum(q0,Q,R);
xhat = x-q0;

%Khat = [system.K,0,0];
%f = system.poly_f(xhat,-Khat*xhat);
K = system.K;
S = system.S;f = system.poly_f(xhat,-K*xhat);
%zero = [0,0;0,0];
%Shat = [system.S,zero;zero,zero];
%V = .5*xhat'*Shat*xhat;
V = .5*xhat'*S*xhat;
%Vdot = diff(V,xhat)*f;
Vdot = diff(V,x)*f;
% found in VDP roa fixed v (SOSTest)
[prog, rho] = prog.newFree(1);
[prog, sig] = prog.newSOSPoly(monomials(x,0:2));
[prog, alph] = prog.newSOSPoly(monomials(x,0:2));
Vdot_sos = (xhat'*xhat)^3*(V - rho) - (1+sig)*Vdot;% - pi*(1-s^2-c^2);
prog = prog.withSOS(Vdot_sos);

% trying to do input saturation
prog = prog.withSOS(-(K*xhat) + 1*alph);

spot_options = spotprog.defaultOptions;
spot_options.verbose = true;
solver = @spot_mosek;
sol = prog.minimize(-rho,solver,spot_options);



N = 300;
[X1,X2] = meshgrid(linspace(-20,20,N), linspace(-20,20,N));
X1DOT = X2;
X2DOT = -system.constants.b*X2-9.81*sin(X1)-(system.K(1).*X1+system.K(2).*X2); 
figure('Position',[0 0 1024 1024])
quiver(X1(:),X2(:),X1DOT(:),X2DOT(:),3.5)
xlim([-20 20])
ylim([-20 20])

hold on
VPLOT = reshape(dmsubs(V,x,[X1(:) X2(:)]'),size(X1));
[~,h] = contour(X1,X2,VPLOT,double(sol.eval(rho))*[1 1]);

set(h,'Color','Red','LineWidth',3)