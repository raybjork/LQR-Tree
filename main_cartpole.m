clc; close all;

%% setup
q0 = [0; 0; 0; 0];
qstar = [0; pi; 0; 0];
Q = diag([10, 10, 1, 1]);
Qf = Q;
R = 1;
system = Cartpole();

N = 10;
u_max = 20;

%% generate trajectory and controller
[x_d, u_d, dt] = collocate_trajectory(system.dynamics(), q0, qstar, u_max, N);
[t, x] = simulate_inputs(system, u, dt, q0);
sytstem.plot(t, x);
[K, S, u] = TVLQR(Q, R, Qf, N * dt, x_d, u_d, u_max, system);

%% simulate and plot
q_err = [0; 0; 0; 0];
f = system.dynamics();
[t, x] = ode45(@(t,x) f(x, u(t,x)), [0 dt*N], x_d(0) + q_err);
system.plot(t, x);