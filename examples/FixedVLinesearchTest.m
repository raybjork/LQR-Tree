%% example from lecture 13
close all
do_continue = true;
rho = .2;
delta = rho*.1;
double_next = true;
rho_opt = 0;
while do_continue
    a = 1;
    prog = spotsosprog();
    [prog,x] = prog.newIndeterminate('x',2);
    f = [-a*(x(1) - x(1)^3/3 - x(2));-x(1)/a];
    
    %% linearization
    A = double(subs(diff(f,x),x,x*0));
    Q = diag([10 1]);
    S = lyap(A',Q);
    
    %% Lyapunov function
    V = .5*x'*S*x;
    Vdot = diff(V,x)*f;
    
    %% SOS constraints
    [prog, sigma_1] = prog.newSOSPoly(monomials(x,6));
    [prog, sigma_2] = prog.newSOSPoly(monomials(x,4));
    Vdot_sos = sigma_1*(V - rho) - (1+sigma_2)*Vdot;
    prog = prog.withSOS(Vdot_sos);
    spot_options = spotprog.defaultOptions;
    spot_options.verbose = true;
    solver = @spot_mosek;
    [prog,slack] = prog.newPos(1);
    sol = prog.minimize(0,solver,spot_options);
    if sol.status ~= spotsolstatus.STATUS_PRIMAL_AND_DUAL_FEASIBLE
        if delta < .05
            do_continue = false;
        else
            delta = delta/2;
            rho = rho - delta;
            double_next = false;
        end
    else
        rho_opt = rho;
        if double_next
            delta = delta*2;
        end
        double_next = true;
        rho = rho + delta;
    end
end

%% PlottingN = 40;
[X1,X2] = meshgrid(linspace(-3,3,N), linspace(-3,3,N));
X1DOT = -a*(X1 - X1.^3/3 - X2);
X2DOT = -X1/a;figure('Position',[0 0 1024 1024])
quiver(X1(:),X2(:),X1DOT(:),X2DOT(:),3.5)
xlim([-3 3])
ylim([-3 3])
x0 = [0;3];
dynamics_forward = @(t,x) [a*(x(1) - x(1)^3/3 - x(2));x(1)/a];
[~,yout] = ode45(dynamics_forward,[0 20],x0);
hold on
plot(yout(50:end,1),yout(50:end,2),'k','LineWidth',3);
VPLOT = reshape(dmsubs(V,x,[X1(:) X2(:)]'),size(X1));
[~,h] = contour(X1,X2,VPLOT,rho_opt*[1 1]);
rho_optset(h,'Color','Red','LineWidth',3)
