function [rho] = TVLQRSOSAlternationsPath(system,AB, S, R, Q, x_d, u_d, dt, N, rho_end)
    rho_init = rho_end;
    state = x_d(0:dt/10:dt*(N-1));
%     rho = ones(N,1)*rho_init*0;
%     return 
    h = plot_end()
    hold on
    plot(state(1,:),state(2,:),'Color','k','LineWidth',3);
    
    hold on
    
    getV = @(x) getLyaps(x,system,AB, S, R,Q, dt,u_d,x_d,N);
    u_max = system.u_max;
    % calculate V and V_dot
    rho_next = rho_end;
    rho = ones(N,1)*rho_init;
    
    
    for i = N:-1:1
        t = (i-1)*dt;
        
        for alt = 1:1
            
            [success,sig,phi,alpha, beta,epsi,gamma,lambda,betsy] = altA(getV,rho(i),rho_next,i,dt,u_max,i == N);
            if ~success
                break;
            end
            [success,rho_out] = altB(getV,sig,phi,alpha,beta,epsi,gamma,lambda,betsy,rho_next,i,dt,u_max,state,i==N,h);
            if ~success
                break;
            end
            rho(i) = double(rho_out);
        end
        rho_next = rho(i);
        if(i>1)
            rho(i-1) = rho(i);
        end
    end
end

function [success,sig,phi,alpha, beta,epsi,gamma,lambda,betsy] = altA(getV,rho,rho_next, i,dt,u_max, last)
    prog = spotsosprog();
    [prog,x] = prog.newIndeterminate('x',2);
    [prog, phi] = prog.newFreePoly(monomials(x,0:4),1);
    [prog, sig] = prog.newSOSPoly(monomials(x,0:4),1);
    [prog, beta] = prog.newSOSPoly(monomials(x,0:4),1);
    [prog, alpha] = prog.newSOSPoly(monomials(x,0:4),1);
    [prog, gamma] = prog.newSOSPoly(monomials(x,0:5),1);
    [prog, lambda] = prog.newSOSPoly(monomials(x,0:5),1);
    [prog, betsy] = prog.newSOSPoly(monomials(x,0:4),1);
    [prog,epsi] = prog.newPos(1);
    [V,V_dot,K_hat] = getV(x);
    
    sol = SOS_constraints(prog,x,V(i),V_dot(i),K_hat(i),rho,rho_next,sig,phi,alpha,beta,epsi,gamma,lambda,betsy,u_max,dt,last);
    phi = sol.eval(phi);
    sig = sol.eval(sig);
    beta = sol.eval(beta);
    alpha = sol.eval(alpha);
    epsi = sol.eval(epsi);
    gamma = sol.eval(gamma);
    lambda = sol.eval(lambda);
    betsy = sol.eval(betsy);
    
    disp(sig)
    
    
    success = sol.status == spotsolstatus.STATUS_PRIMAL_AND_DUAL_FEASIBLE;
end

function [success,rho] = altB(getV,sig,phi,alpha,beta,epsi,gamma,lambda,betsy,rho_next,i,dt,u_max,state,last,h)
    prog = spotsosprog();
    [prog,x] = prog.newIndeterminate('x',2);
    [prog, rho] = prog.newPos(1);
    [V,V_dot,K_hat] = getV(x);
    
    sol = SOS_constraints(prog,x,V(i),V_dot(i),K_hat(i),rho,rho_next,sig,phi,alpha,beta,epsi,gamma,lambda,betsy,u_max,dt,last);
    rho = sol.eval(rho);
    
    success = sol.status == spotsolstatus.STATUS_PRIMAL_AND_DUAL_FEASIBLE;
    if success
        PlotFunnel(rho,V(i),x,state,h);
    end
end

function sol = SOS_constraints(prog,x,V,V_dot,K_hat,rho,rho_next,sig,phi,alpha,beta,epsi,gamma,lambda,betsy,u_max,dt,last)
    
    if(isa(rho,'double'))
        prog = prog.withSOS((beta)*(V - rho) - (gamma)*(-K_hat-u_max)- epsi * (x' * x)^2);
        prog = prog.withSOS((betsy)*(V - rho) - (lambda)*(K_hat-u_max)- epsi * (x' * x)^2);
    else
        prog = prog.withSOS((beta)*(V - rho) - (gamma)*(-K_hat-u_max));
        prog = prog.withSOS((betsy)*(V - rho) - (lambda)*(K_hat-u_max));
    end
    %prog = prog.withSOS(V);
    
    
        %Forward derivative
        rho_dot = (rho_next-rho)/dt;
        %ensure V_dot is less than rho_dot when v=rho
        if(isa(rho,'double'))
            prog = prog.withSOS((sig)*(rho_dot - V_dot) + phi*(V-rho) - epsi * (x' * x)^2);
            %prog = prog.withSOS((sig)*(rho_dot - V_dot) + phi*(V-rho));
                
                        
            %prog = prog.withSOS((phi(i))*(V(i)-rho(i)));
            %prog = prog.withSOS((1+sig(i)) - V_dot(i) + rho_dot + phi(i)*(V(i)-rho(i)));
        else
            prog = prog.withSOS((sig)*(rho_dot - V_dot) + phi*(V-rho));
            %prog = prog.withSOS((sig)*(rho_dot - V_dot) + phi*(V-rho)); 
                    
                    
        end
        
    if(last)
        if(isa(rho,'double'))
            prog = prog.withSOS(alpha*(rho_next - rho));
        else
            prog = prog.withSOS(alpha*(rho_next - rho));
        end
    end
    spot_options = spotprog.defaultOptions;
    spot_options.verbose = true;
    solver = @spot_mosek;
    if(isa(rho,'double'))
        sol = prog.minimize(-epsi,solver,spot_options);
    else
        sol = prog.minimize(-rho,solver,spot_options);
    end
    if sol.status ~= spotsolstatus.STATUS_PRIMAL_AND_DUAL_FEASIBLE
        
        
%         for r = 0:.5:2
%             figure
%             title(['rho = ', num2str(r)]);
%             lowest = inf;
%             for x2 = 0:.05:6
%                 for x1 = -pi:.05:6
% 
%                     hold on
%                     val = subs((sig)*(rho_dot - V_dot) + phi*(V-rho),x(1),x1);
%                     val = subs(val,x(2),x2);
%                     val = subs(val,rho,r);
%                     if(double(val) < lowest)
%                         lowest = double(val);
%                     end
%                     %plot3(x1,x2,val,'*');
%                 end
%             end
%             disp(lowest)
%          end
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

function [] = PlotFunnel(rho, V,x,state,h)
    

    sos_color = 1/255*[131 223 255];
    plot_size = 6;
    density = 300;
    [X1,X2] = meshgrid(linspace(-plot_size,plot_size,density), linspace(-plot_size,plot_size,density));
    VPLOT = reshape(dmsubs(-V,x,[X1(:) X2(:)]'),size(X1));
    [~,sos] = contourf(X1,X2,VPLOT,-double(rho)*[1 1]);
    colormap(sos_color)

    set(sos,'Color',sos_color*.9,'LineWidth',3)
    
    
    traj = plot(state(1,:),state(2,:),'Color','b','LineWidth',2);
    hold on
   
    
    l = legend([sos,traj],'Infinite-Time LQR ROA','SOS ROA','Nominal Trajectory');
    set(l,'Location','southwest');
    
    title('SOS ROA For Nominal Trajectory: Individual Alternations','interpreter','latex');
    xlabel('$\theta$','interpreter','latex');
    ylabel('$\dot{\theta}$','interpreter','latex');
    
end

function [h] = plot_end()
    
     load infinite_SOS_data.mat
    N = 300;
    [X1,X2] = meshgrid(linspace(0,6,N), linspace(-6,6,N));

    VPLOT = reshape(dmsubs(-V,x,[X1(:) X2(:)]'),size(X1));
    [~,h] = contour(X1,X2,VPLOT,-double(rho)*[1 1]);
    set(h,'Color','k','LineWidth',3)
end