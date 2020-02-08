clc; close all; clear
addpath(genpath(pwd))
%% setup
qstar = [pi;0]; % fixed point for system

Q = diag([1 1]); % initial Q
R = 1; % initial R

u_max = 5;

system = Pendulum(qstar, Q, R, u_max); % infinite time LQR controller

N = 50; % number of collocation points

M = 10; % number of nodes to sample in tree 
q_max = [2*pi; 3]; % limits of LQR tree exploration


%% build tree
r = infinite_SOS(system);
paths(1) = Trajectory(system.S, @(t)r, 1, 1, qstar, 0); % infinite time to seed tree
qset = [-2,0
        0,0; 3.5, -4]
%qset = [0,0];      
        
for i = 2:M+1
    point_valid = false;
%     while ~point_valid
%         q = [rand * q_max(1); rand * q_max(2)];
%         point_valid = check_startpoint(q,paths);
        q = qset(i-1,:)';
%     end
    closest_path = Inf;
    for j = 1:i-1
        for t = 0 :paths(j).dt: paths(j).dt*(paths(j).N-1)
            %V = (q-paths(j).x_d(t))' * paths(j).S(t) * (q-paths(j).x_d(t));
            V = Tedrake_proximity(q, paths(j).x_d(t),R, system,3);
            if V < closest_path(1)
                closest_path = [V j t];
            end
        end
    end
    rho_end = paths(closest_path(2)).rho(closest_path(3));
    if (q-paths(j).x_d(t))' * paths(closest_path(2)).S(closest_path(3)) * (q-paths(j).x_d(t)) >= rho_end
        % generate trajectory and controller
        [x_d, u_d, dt] = collocate_trajectory(q, paths(closest_path(2)).x_d(closest_path(3)), N, system);
        
        
        %save(['LQRTreeTraj',num2str(i),'.mat'] , 'x_d', 'u_d', 'dt', 'N','system');
        [S, AB, u] = TVLQR(x_d, u_d, dt * (N-1),paths(closest_path(2)).S(closest_path(3)), system);

        % simulate and plot
        q_err = [0; 0];
        f = system.dynamics();
        [t, x] = ode45(@(t,x) f(x, u(t,x)), [0 dt*(N-1)], x_d(0) + q_err);
        %system.plot(t, x);
        %%
        % alternations to find ROA of trajectory
       
        rho = TVLQRSOSAlternationsPath(system, AB, S, R, Q, x_d, u_d, dt, N, rho_end);
        
        rho = @(t) ppval(spline(0:dt:dt*(N-1), rho), t);
        %rho = TVLQRSOS(system,AB, S, R, Q, u_max, x_d, u_d, dt, N);
        
        % append new trajectory
        state = x_d(0:dt:(N-1)*dt);
        plot(state(1,:), state(2,:));
        paths(i) = Trajectory(S, rho, N, dt, x_d, u_d);
        drawnow;
    end
end

function [valid] = check_startpoint(q,paths)
    valid = true;
    for i = 1:size(paths,2)
        dt = paths(i).dt;
        N = paths(i).N;
        x_d = paths(i).x_d;
        S = paths(i).S;
        rho = paths(i).rho;
        
        for t = 0:dt:dt*(N-1)
            test_rho = (q-x_d(t))'*S(t)*(q-x_d(t));
            if test_rho < rho(t)
                valid = false;
                break;
            end
        end
        
        if ~valid
            break;
        end
    end
   
end