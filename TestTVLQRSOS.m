clc; close all;

%% setup
q0 = [0; 0];
qstar = [pi;0];
Q = diag([1 1]);
Qf = diag([1 1]);
R = 1;
system = Pendulum(qstar,Q,R);

N = 40;
u_max = 20;

%% generate trajectory and controller
[x_d, u_d, dt] = collocate_trajectory(system.dynamics(), q0, qstar, u_max, N);
[S ,AB, u] = TVLQR(Q, R, Qf, N * dt, x_d, u_d,u_max, system);

%% plot SOS funnel
TVLQRSOSAlternations(system, AB,S,R,Q, u, u_max,x_d,u_d,dt,N);