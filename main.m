clc; close all;

%% setup
q0 = [0; 0];
qstar = [pi;0];
Q = diag([1 1]);
Qf = diag([1 1]);
R = 1;
system = Pendulum(qstar,Q,R);

N = 10;
u_max = 20;

%% generate trajectory and controller
[x_d, u_d, u_t, dt] = collocate_trajectory(system.dynamics(), q0, qstar, u_max, N);
[K, u] = TVLQR(Q, R, Qf, N * dt, x_d, u_d,u_max, system);

%% simulate and plot
q_err = [2; 0];
f = system.dynamics();
[t, x] = ode45(@(t,x) f(x, u(t,x)), [0 dt*N], ppval(x_d, 0) + q_err);
system.plot(t, x);