function x = TVLQRSafe(system,AB, S, R, Q, u, u_max, x_d, u_d, dt, N)

state = x_d(0:dt/10:dt*N);
figure
hold on
plot(state(1,:),state(2,:),'Color','k','LineWidth',3);

hold on

for t = 0:dt:dt*N
    q_0 = x_d(t);
    prog = spotsosprog();
    [A,B] = AB(t);
    K = inv(R)*B'*S(t);
    [prog,x] = prog.newIndeterminate('x',2);
    
    x_hat = x-q_0;
    f = system.poly_f(x, u_d(t)-K*x_hat);
    V = .5*x_hat'*S(t)*x_hat;
    
    %need to add ricatti equation because chain rule
    S_dot = -S(t)*A +A'*S(t)+Q-S(t)*B*inv(R)*B'*S(t);%Ricatti
    V_dot = diff(V,x)*f + x_hat'*S_dot*x_hat;
    
    [prog, rho] = prog.newFree(1);
    [prog, sig] = prog.newSOSPoly(monomials(x,0:2));
    [prog, alph] = prog.newSOSPoly(monomials(x,0:2));
    [prog, beta] = prog.newSOSPoly(monomials(x,0:2));
    
    Vdot_sos = (x_hat'*x_hat)^3*(V - rho) - (1+sig)*V_dot ;% - pi*(1-s^2-c^2);
   
    prog = prog.withSOS(Vdot_sos);

    %prog = prog.withSOS((V - rho) - (1+beta)*(-K(t)*xhat-u_max));
    %prog = prog.withSOS((V - rho) - (1+alph)*(K(t)*xhat-u_max));


    spot_options = spotprog.defaultOptions;
    spot_options.verbose = true;
    solver = @spot_mosek;
    sol = prog.minimize(-rho,solver,spot_options);
    
    
    plot_size = 6;
    density = 300;
    [X1,X2] = meshgrid(linspace(-plot_size,plot_size,density), linspace(-plot_size,plot_size,density));
    VPLOT = reshape(dmsubs(V,x,[X1(:) X2(:)]'),size(X1));
    [~,h] = contour(X1,X2,VPLOT,double(sol.eval(rho))*[1 1]);
    set(h,'Color','Red','LineWidth',3)
end


for i = -1:1:1
    for j = -1:1:1
        init_state = x_d(0);
        init_state = [init_state(1)+i; init_state(2)+j];
        f = system.dynamics();
        [~, x] = ode45(@(t,x) f(x, u(t,x)), [0 dt*N], init_state);
        plot(x(:,1),x(:,2),'Color','g');
    end
end


