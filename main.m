close all; clear all; clc;

q0 = [0; 0];
qstar = [pi;0];
Q = diag([1 1]);
Qf = diag([0 0]);
R = 1;
system = Pendulum(qstar,Q,R);

N = 30;
dt = 6/N;
u_max = 20;

[x, u] = collocate_trajectory(system.dynamics(), q0, qstar, u_max, N, dt);

systems = system.empty(length(x), 0);
for i = 1:length(x)
    systems(i) = Pendulum(x(i,:)', Q, R);
end
opts = odeset('MaxStep', 0.01);
[t, S] = ode45(@(t,S) TV_Ricatti(t, S, systems), [0 dt*N], Qf, opts);


[t, x] = simulate_inputs(system, u, dt, q0);
system.plot(t, x);


