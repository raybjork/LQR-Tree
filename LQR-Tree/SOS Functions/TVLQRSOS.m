function rho = TVLQRSOS(system,AB, S, R, Q, u_max, x_d, u_d, dt, N)
    
    state = x_d(0:dt/10:dt*N);
    figure
    hold on
    plot(state(1,:),state(2,:),'Color','k','LineWidth',3);
    
    hold on
    prog = spotsosprog();
    [prog,x] = prog.newIndeterminate('x',2);
    [prog, rho] = prog.newFree(1);
    V = zeros(1,N)*x(1);
    
    for t = 0:dt:dt*N
        i = t/dt+1;
        
        q_0 = x_d(t);
        [A,B] = AB(t);
        K = inv(R)*B'*S(t);
        
        
        x_hat = x-q_0;
        f = system.poly_f(x, u_d(t)-K*x_hat,x_d(t));
        V(i) = .5*x_hat'*S(t)*x_hat;
        
        %need to add ricatti equation because chain rule
        S_dot = -(S(t)*A +A'*S(t)+Q-S(t)*B*inv(R)*B'*S(t));%Ricatti
        V_dot = diff(V(i),x)*f + .5*x_hat'*S_dot*x_hat;
        
        %[prog, beta] = prog.newSOSPoly(monomials(x,0:2));
        
        if(i<N)
            [prog, sig] = prog.newSOSPoly(monomials(x,0:4));
            %Forward derivative
            
            
            %ensure V_dot is less than rho_dot when v=rho
            %prog = prog.withSOS((1+sig)*(rho_dot - V_dot) + phi*(V-rho(i)));
            Vdot_sos = (x_hat'*x_hat)^3*(V(i) - rho) - (1+sig)*V_dot;
            prog = prog.withSOS(Vdot_sos);
        else
            %[prog, beta] = prog.newFreePoly(monomials(x,0:4),1);
            %         [prog, alpha] = prog.newSOSPoly(monomials(x,0:4),1);
            %         prog = prog.withSOS(alpha*(rho - V) - (rho-1000));
            %
        end
        
        
        %prog = prog.withSOS((V - rho) - (1+beta)*(-K(t)*xhat-u_max));
        %prog = prog.withSOS((V - rho) - (1+alph)*(K(t)*xhat-u_max));
        
        
        
    end
    
    spot_options = spotprog.defaultOptions;
    spot_options.verbose = true;
    solver = @spot_mosek;
    sol = prog.minimize(-rho,solver,spot_options);
    
    rho = sol.eval(rho);
    PlotFunnel(rho,V,x,N);
    
    
    %     plot_size = 6;
    %     density = 300;
    %     [X1,X2] = meshgrid(linspace(-plot_size,plot_size,density), linspace(-plot_size,plot_size,density));
    %     VPLOT = reshape(dmsubs(V,x,[X1(:) X2(:)]'),size(X1));
    %     [~,h] = contour(X1,X2,VPLOT,double(sol.eval(rho))*[1 1]);
    %     set(h,'Color','Red','LineWidth',3)
    
end

function [] = PlotFunnel(rho, V,x, N)
    
    hold on
    for i = 1:N
        plot_size = 6;
        density = 300;
        [X1,X2] = meshgrid(linspace(-plot_size,plot_size,density), linspace(-plot_size,plot_size,density));
        VPLOT = reshape(dmsubs(V(i),x,[X1(:) X2(:)]'),size(X1));
        [~,h] = contour(X1,X2,VPLOT,double(rho)*[1 1]);
        set(h,'Color','Red','LineWidth',3)
        
    end
    
end