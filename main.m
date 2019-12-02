clc; close all;

%% setup
q0 = [0; 0];
qstar = [pi;0];
Q = diag([1 1]);
Qf = diag([1 1]);
R = 1;

N = 10;
u_max = 20;

system = Pendulum(qstar,Q,R, u_max);

%% generate trajectory and controller
[x_d, u_d, dt] = collocate_trajectory(q0, qstar, N, system);
[K, S, u] = TVLQR(x_d, u_d, dt * N, system);

%% simulate and plot
q_err = [1; 0];
f = system.dynamics();
[t, x] = ode45(@(t,x) f(x, u(t,x)), [0 dt*N], x_d(0) + q_err);
system.plot(t, x);