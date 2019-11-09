system = Cartpole();

Q = diag([1 1 1 10]);
R = 1;

[K, S] = lqr(system, Q, R);

tf = 10;
x_0 = [0; pi+1; 0; 0];

[t, x] = system.simulate(K, tf, x_0);
system.plot(t, x)
