clc; close all;

%% setup
q0 = [0; 0];        % initial conditions
qstar = [pi;0];     % fixed point
Q = diag([1 1]);
R = 1;             
u_max = 20;         % torque limit

system = Pendulum(qstar, Q, R, u_max);

N = 10;     % number of collocation points

%% generate trajectory and controller
[x_d, u_d, dt] = collocate_trajectory(qstar, q0, N, system);
[S, AB, u] = TVLQR(x_d, u_d, N * dt, system);

%% simulate and plot
q_err = [1; -3];
f = system.dynamics();
[t, x] = ode45(@(t,x) f(x, u(t,x)), [0 dt*N*1.5], x_d(0) + q_err);
%system.plot(t, x);
syms x [2 1] real
getLyapstrue(x, system, AB, S,R, Q, dt, u_d, x_d, N);