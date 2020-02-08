function [rho] = TVLQRSOSAlternations(system,AB, S, R, Q, x_d, u_d, dt, N, rho_end)
    rho_init = .2;
    state = x_d(0:dt/10:dt*(N-1));
    
    hold on
    plot(state(1,:),state(2,:),'Color','k','LineWidth',3);
    
    hold on
    rho = ones(1,N)*rho_init;
    getV = @(x) getLyaps(x,system,AB, S, R,Q, dt,u_d,x_d,N);
    u_max = system.u_max;
    % calculate V and V_dot
    for i = 1:5
        
        [success,sig,phi,alpha, beta,epsi,gamma,lambda] = altA(getV,rho,dt,N,u_max,rho_end);
        if ~success
            break;
        end
        [success,rho] = altB(getV,sig,phi,alpha,beta,epsi,gamma,lambda,dt,N,u_max,state,rho_end);
        if ~success
            break;
        end
    end   
    rho = double(rho);
end

function [success,sig,phi,alpha, beta,epsi,gamma,lambda] = altA(getV,rho,dt,N,u_max, rho_end)
    prog = spotsosprog();
    [prog,x] = prog.newIndeterminate('x',2);
    [prog, phi] = prog.newFreePoly(monomials(x,0:4),N);
    [prog, sig] = prog.newSOSPoly(monomials(x,0:4),N);
    [prog, beta] = prog.newFreePoly(monomials(x,0:4),1);
    [prog, alpha] = prog.newSOSPoly(monomials(x,0:4),1);
    [prog, gamma] = prog.newSOSPoly(monomials(x,0:4),N);
    [prog, lambda] = prog.newSOSPoly(monomials(x,0:4),N);
    [prog,epsi] = prog.newPos(1);
    [V,V_dot,K_hat] = getV(x);
    sol = SOS_constraints(prog,x,V,V_dot,K_hat,rho,sig,phi,alpha,beta,epsi,gamma,lambda,dt,N,u_max,rho_end);
    phi = sol.eval(phi);
    sig = sol.eval(sig);
    beta = sol.eval(beta);
    alpha = sol.eval(alpha);
    epsi = sol.eval(epsi);
    gamma = sol.eval(gamma);
    lambda = sol.eval(lambda);


    success = sol.status == spotsolstatus.STATUS_PRIMAL_AND_DUAL_FEASIBLE;
end

function [success,rho] = altB(getV,sig,phi,alpha,beta,epsi,gamma,lambda,dt,N,u_max,state,rho_end)
    prog = spotsosprog();
    [prog,x] = prog.newIndeterminate('x',2);
    [prog, rho] = prog.newPos(N);
    [V,V_dot,K_hat] = getV(x);
    sol = SOS_constraints(prog,x,V,V_dot,K_hat,rho,sig,phi,alpha,beta,epsi,gamma,lambda,dt,N,u_max,rho_end);
    rho = sol.eval(rho);
    
    success = sol.status == spotsolstatus.STATUS_PRIMAL_AND_DUAL_FEASIBLE;
    if success
        PlotFunnel(rho,V,x,N,state);
    end
end

function sol = SOS_constraints(prog,x,V,V_dot,K_hat,rho,sig,phi,alpha,beta,epsi,gamma,lambda,dt,N,u_max,rho_end)
    for i = 1:N
        
        prog = prog.withSOS((V(i) - rho(i)) - (1+gamma(i))*(-K_hat(i)-u_max));
        prog = prog.withSOS((V(i) - rho(i)) - (1+lambda(i))*(K_hat(i)-u_max));
        if (i==1)
            %prog=prog.withSOS(rho(i)-V(i));
        end
       prog = prog.withSOS(V(i));
        
        if(i<N)
            %Forward derivative
            rho_dot = (rho(i+1)-rho(i))/dt;
            %ensure V_dot is less than rho_dot when v=rho
            if(isa(rho,'double'))
                prog = prog.withSOS((sig(i))*(rho_dot - V_dot(i)) + phi(i)*(V(i)-rho(i)) - epsi * (x' * x)^2); 
                %prog = prog.withSOS((phi(i))*(V(i)-rho(i)));
                %prog = prog.withSOS((1+sig(i)) - V_dot(i) + rho_dot + phi(i)*(V(i)-rho(i)));
            else
                prog = prog.withSOS((sig(i))*(rho_dot - V_dot(i)) + phi(i)*(V(i)-rho(i)) -epsi * (x' * x)^2); 
            end
            
        else
            if(isa(rho,'double'))
                prog = prog.withSOS(-alpha*(rho(i) - V(i)) - beta*(rho(i)-rho_end));
            else
                prog = prog.withSOS(-alpha*(rho(i) - V(i)) - beta*(rho(i)-rho_end));
            end
        end
    end 
    spot_options = spotprog.defaultOptions;
    spot_options.verbose = true;
    solver = @spot_mosek;
    if(isa(rho,'double'))
        sol = prog.minimize(-epsi,solver,spot_options);
    else
        sol = prog.minimize(-rho(N),solver,spot_options);
    end
end

function [V,V_dot,K_hat] = getLyaps(x,system, AB, S, R,Q, dt,u_d,x_d,N)
    V = zeros(1,N)*x(1);
    V_dot = zeros(1,N)*x(1);
    K_hat = zeros(1,N)*x(1);
    
    for i = 1:N
        t = (i-1)*dt;
        q_0 = x_d(t);
        [A,B] = AB(t);
        K = R\B'*S(t);

        x_hat = x-q_0;
        K_hat(i) = K*x_hat;
        %f = system.poly_f(x, u_d(t)-K*x_hat);
        %this is significantly better because V_dot is negative definite
        %everywhere 
        f = system.poly_f(x,u_d(t)-K*x_hat,q_0);
        V(i) = .5*x_hat'*S(t)*x_hat;
        
        %need to add ricatti equation because chain rule
        S_dot = -(S(t)*A +A'*S(t)+Q-S(t)*B*inv(R)*B'*S(t));%Ricatti
        V_dot(i) = diff(V(i),x)*f + .5*x_hat'*S_dot*x_hat;
        %v = eig(double(subs(diff(diff(V_dot(i), x)', x), x, q_0)))
        
    end
    
    
end

function [] = PlotFunnel(rho, V,x, N,state)
    save trajectorySOSdata.mat rho V x N state
    for i = 1:N
        plot_size = 6;
        density = 300;
        [X1,X2] = meshgrid(linspace(-plot_size,plot_size,density), linspace(-plot_size,plot_size,density));
        VPLOT = reshape(dmsubs(-V(i),x,[X1(:) X2(:)]'),size(X1));
        [~,h] = contourf(X1,X2,VPLOT,-double(rho(i))*[1 1]);
        set(h,'Color','Red','LineWidth',3)
        
    end
    
    plot(state(1,:),state(2,:),'Color','k','LineWidth',3);
    hold on

end