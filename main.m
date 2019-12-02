clc; clf;
drawnow;
hold on;

%% setup
qstar = [pi;0]; % fixed point for system

Q = diag([1 1]); % initial Q
R = 1; % initial R

u_max = 10000;

system = Pendulum(qstar, Q, R, u_max); % infinite time LQR controller

N = 40; % number of collocation points

M = 30; % number of nodes to sample in tree 
%paths = zeros(M + 1, 1); % list of nodes for tree
q_max = [2*pi; 3]; % limits of LQR tree exploration

%% build tree
%rho = infinite_SOS(system)
rho = 0;
paths(1) = Trajectory(system.S, rho, 1, 1, qstar); % infinite time to seed tree

A = zeros(N + 1);
A(1,1) = 0;

for i = 2:M+1
    q = [rand * q_max(1); rand * q_max(2)];
    sample = Pendulum(q, Q, R, u_max);
    closest_path = Inf;
    for j = 1:i-1
        for t = 0 :paths(j).dt: paths(j).dt*N
            V = (q-paths(j).x_d(t))' * sample.S * (q-paths(j).x_d(t));
            if V < closest_path(1)
                closest_path = [V j t];
            end
        end
    end
    q_connection = paths(closest_path(2)).x_d(closest_path(3));
    if q' * paths(closest_path(2)).S(closest_path(3)) * q >= 0 %root.rho
        % generate trajectory and controller
        [x_d, u_d, dt] = Copy_of_collocate_trajectory(q, q_connection, N, system);
        [S, AB, u] = TVLQR(x_d, u_d, dt * N,paths(closest_path(2)).S(closest_path(3)), system);
        
        % simulate and plot
        q_err = [0; 0];
        f = system.dynamics();
        [t, x] = ode45(@(t,x) f(x, u(t,x)), [0 dt*N], x_d(0) + q_err);
        %system.plot(t, x);
        
        
        % alternations to find ROA of trajectory
        %rho = TVLQRSOSAlternations(system, AB, S, R, Q, x_d, u_d, dt, N, root.rho);

        % append new trajectory
        paths(i) = Trajectory(S, rho*0, N, dt, x_d);
        
%         cost = (q - q_connection)' * S(0) * (q - q_connection);
%         A(i, closest_path(2)) = cost;
%         A(closest_path(2), i) = cost;
%         A(i, i) = A(closest_path(2), closest_path(2)) + cost;
        
        t = 0 : dt : dt * (N-1);
        states = paths(i).x_d(t);
        plot(states(1, :), states(2, :), 'LineWidth', 3)
        drawnow
        
    end
end