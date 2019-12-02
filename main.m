clc; close all;

%% setup
qstar = [pi;0]; % fixed point for system

Q = diag([1 1]); % initial Q
R = 1; % initial R

u_max = 50;

system = Pendulum(qstar, Q, R, u_max); % infinite time LQR controller

N = 40; % number of collocation points

M = 10; % number of nodes to sample in tree 
%paths = zeros(M + 1, 1); % list of nodes for tree
q_max = [2*pi; 3]; % limits of LQR tree exploration


%% build tree
paths(1) = Trajectory(system.S, infinite_SOS(system), 1, 1, qstar); % infinite time to seed tree
for i = 2:M+1
    q = [rand * q_max(1); rand * q_max(2)];
    closest_path = Inf;
    for j = 1:i-1
        for t = 0 :dt: dt*N
            V = (q-paths(j).x_d(t))' * paths(j).S(t) * (q-paths(j).x_d(t));
            if V < closest_path(1)
                closest_path = [V j t];
            end
        end
    end
    closest_path
    if q' * paths(closest_path(2)).S(closest_path(3)) * q >= 0 %root.rho
        % generate trajectory and controller
        [x_d, u_d, dt] = collocate_trajectory(q, paths(closest_path(2)).x_d(closest_path(3)), N, system);
        [S, AB, u] = TVLQR(x_d, u_d, dt * N,paths(closest_path(2)).S(closest_path(3)), system);
        
        % simulate and plot
        q_err = [0; 0];
        f = system.dynamics();
        [t, x] = ode45(@(t,x) f(x, u(t,x)), [0 dt*N], x_d(0) + q_err);
        %system.plot(t, x);
        
        
        % alternations to find ROA of trajectory
        rho = TVLQRSOSAlternations(system, AB, S, R, Q, x_d, u_d, dt, N, root.rho);

        % append new trajectory
        paths(i) = Trajectory(S, rho*0, N, dt, x_d);
        
    end
end