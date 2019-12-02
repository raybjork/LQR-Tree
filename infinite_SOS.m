function [rho] = infinite_SOS(system)
prog = spotsosprog();

[prog,x] = prog.newIndeterminate('x',2);


xhat = x-system.qstar;

K = inv(system.R)*system.B'*system.S;
S = system.S;
f = system.poly_f(xhat,-K*xhat,system.qstar);
V = .5*xhat'*S*xhat;
Vdot = diff(V,x)*f;

% found in VDP roa fixed v (SOSTest)
[prog, rho] = prog.newFree(1);
[prog, sig] = prog.newSOSPoly(monomials(x,0:4));
[prog, alpha] = prog.newSOSPoly(monomials(x,0:4));
[prog, beta] = prog.newSOSPoly(monomials(x,0:4));
Vdot_sos = (xhat'*xhat)^3*(V - rho) - (1+sig)*Vdot;
prog = prog.withSOS(Vdot_sos);

% trying to do input saturation
prog = prog.withSOS((V - rho) - (1+beta)*(-K*xhat-system.u_max));
prog = prog.withSOS((V - rho) - (1+alpha)*(K*xhat-system.u_max));

spot_options = spotprog.defaultOptions;
spot_options.verbose = true;
solver = @spot_mosek;
sol = prog.minimize(-rho,solver,spot_options);
rho = double(sol.eval(rho));

N = 300;
[X1,X2] = meshgrid(linspace(-2*pi,2*pi,N), linspace(-2*pi,2*pi,N));


hold on
VPLOT = reshape(dmsubs(V,x,[X1(:) X2(:)]'),size(X1));
[~,h] = contour(X1,X2,VPLOT,double(sol.eval(rho))*[1 1]);

set(h,'Color','Red','LineWidth',3)
end

