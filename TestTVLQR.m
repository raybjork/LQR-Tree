q0 = [0; 0];
qstar = [pi;0];
Q = diag([1 1]);
R = 1;
system = Pendulum(qstar,Q,R);

N = 30;
dt = 6/N;
u_max = 20;

[xd,ud] = collocate_trajectory(system.dynamics(), q0, qstar, u_max, N, dt);
[t, x] = SimulateTVLQR(system, u, dt, q0);
system.plot(t, x);