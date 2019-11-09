%% example from lecture 13
close all 
clear 
a = 1;
prog = spotsosprog();

[prog,x] = prog.newIndeterminate('x',2);
f = [-a*(x(1) - x(1)^3/3 - x(2)); -x(1)/a];
pt = [1;1];
%% linearization
A = double(subs(diff(f,x),x,x*0));
Q = diag([5 1]);
S = lyap(A',Q);

%% Lyapunov function
V = .5*x'*S*x;
Vdot = diff(V,x)*f;
%% SOS constraints
[prog, rho] = prog.newFree(1);
[prog, sigma_2] = prog.newSOSPoly(monomials(x,4));
Vdot_sos = (x'*x)^3*(V - rho) - (1+sigma_2)*Vdot;
prog = prog.withSOS(Vdot_sos);
spot_options = spotprog.defaultOptions;
spot_options.verbose = true;
solver = @spot_mosek;
sol = prog.minimize(-rho,solver,spot_options);

%% Plotting
N = 20;
[X1,X2] = meshgrid(linspace(-3,3,N), linspace(-3,3,N));
X1DOT = -a*(X1 - X1.^3/3 - X2);
X2DOT = -X1/a;
figure('Position',[0 0 1024 1024])
quiver(X1(:),X2(:),X1DOT(:),X2DOT(:),3.5)
xlim([-3 3])
ylim([-3 3])
x0 = [0;3];
dynamics_forward = @(t,x) [a*(x(1) - x(1)^3/3 - x(2)); x(1)/a];
[~,yout] = ode45(dynamics_forward,[0 10],x0);
hold on
plot(yout(50:end,1),yout(50:end,2),'k','LineWidth',3);
VPLOT = reshape(dmsubs(V,x,[X1(:) X2(:)]'),size(X1));
[~,h] = contour(X1,X2,VPLOT,double(sol.eval(rho))*[1 1]);
set(h,'Color','Red','LineWidth',3)
