%% example from Lecture 13

function vdp_roa_alternations
close all

degree = 2;

a = 1;
prog = spotsosprog();
[prog,x] = prog.newIndeterminate('x',2);

f = [-a*(x(1) - x(1)^3/3 - x(2)); -x(1)/a];

%% linearization
A = double(subs(diff(f,x),x,x*0));
Q = diag([1 1]);
S = lyap(A',Q);

%% Lyapunov function
V = .5*x'*S*x;
rho = .01;

V = V/double(subs(V,x,[1;1]));


plotVectorField

for iter=1:20
  [success,sigma_1, sigma_2] = alternationA(V,rho);
  if ~success
    do_continue = false;
    break;
  end
  [success,V, rho] = alternationB(sigma_1,sigma_2, degree);
  if ~success
    do_continue = false;
    break;
  end
  
  plotV(V,x,rho);
%   pause
end
end

%% Alternation A
function [success,sigma_1, sigma_2] = alternationA(V,rho)
a = 1;
prog = spotsosprog();
[prog,x] = prog.newIndeterminate('x',2);

f = [-a*(x(1) - x(1)^3/3 - x(2)); -x(1)/a];

Vdot = diff(V,x)*f;


%% SOS constraints
% define sigma_1 and sigma_2 to be homoegenous degree 4 and 2 polynomials
[prog, sigma_1] = prog.newSOSPoly(monomials(x,2:4));
[prog, sigma_2] = prog.newSOSPoly(monomials(x,0:2));

[prog, gamma] = prog.newPos(1);

Vdot_sos = sigma_1*(V - rho) - (1+sigma_2)*Vdot;
prog = prog.withSOS(Vdot_sos - gamma*x'*x);

spot_options = spotprog.defaultOptions;
spot_options.verbose = true;
solver = @spot_mosek;
sol = prog.minimize(-gamma,solver,spot_options);

% back off--solve a feasability problem (no objective)
% where gamma > 0.98 gamma* (gamma* being the optimal value)
% this is a numerical "trick" to help improve the alternations

%temporarily removed
 prog = prog.withPos(gamma - .98*sol.eval(gamma));
 
 sol = prog.minimize(0*x(1),solver,spot_options);

sigma_1 = sol.eval(sigma_1);
sigma_2 = sol.eval(sigma_2);

success = sol.status == spotsolstatus.STATUS_PRIMAL_AND_DUAL_FEASIBLE;
end

%% Alternation B
function [success,V, rho] = alternationB(sigma_1,sigma_2, degree)
a = 1;
prog = spotsosprog();
[prog,x] = prog.newIndeterminate('x',2);

f = [-a*(x(1) - x(1)^3/3 - x(2)); -x(1)/a];

%% Lyapunov function
% define V to be a homogeneous polynomial of given degree. A more general
% choice would be to use monomials(x,0:degree) which would include all
% monomials of UP TO the given degree
[prog, V] = prog.newSOSPoly(monomials(x,2:degree));
rho = 1;
% [prog, rho] = prog.newFree(1);
% prog = prog.withEqs(subs(V,x,[1;1])-1); % Fix V(1,1) = 1;
% prog = prog.withEqs(subs(V,x,[-1;-1])-1); % Fix V(1,1) = 1;
cost = subs(V,x,[1;1]);% + subs(V,x,[-1;-1]);% + subs(V,x,[-1;1]) + subs(V,x,[1;-1]);
Vdot = diff(V,x)*f;


%% SOS constraints

Vdot_sos = sigma_1*(V - rho) - (1+sigma_2)*Vdot;
prog = prog.withSOS(Vdot_sos);

spot_options = spotprog.defaultOptions;
spot_options.verbose = true;
solver = @spot_mosek;
sol = prog.minimize(cost,solver,spot_options);
success = sol.status == spotsolstatus.STATUS_PRIMAL_AND_DUAL_FEASIBLE;
if ~success
  return;
end
V = sol.eval(V);
rho = double(sol.eval(rho));

% rescale the solution to use rho=1
V = V/rho;
rho = 1;


end

%% Plotting
function plotVectorField
a = 1;
N = 60;
[X1,X2] = meshgrid(linspace(-3,3,N), linspace(-3,3,N));
X1DOT = -a*(X1 - X1.^3/3 - X2);
X2DOT = -X1/a;

h = figure(1);
set(h,'Position',[0 0 1024 1024])
quiver(X1(:),X2(:),X1DOT(:),X2DOT(:),3.5)

xlim([-3 3])
ylim([-3 3])

x0 = [0;3];
dynamics_forward = @(t,x) [a*(x(1) - x(1)^3/3 - x(2)); x(1)/a];
[~,yout] = ode45(dynamics_forward,[0 20],x0);
hold on
plot(yout(50:end,1),yout(50:end,2),'k','LineWidth',3);

end

function plotV(V,x,rho)
N = 30;
[X1,X2] = meshgrid(linspace(-3,3,N), linspace(-3,3,N));
VPLOT = reshape(dmsubs(V,x,[X1(:) X2(:)]'),size(X1));

[~,h] = contour(X1,X2,VPLOT,rho*[1 1]);

set(h,'Color','Red','LineWidth',3)
end