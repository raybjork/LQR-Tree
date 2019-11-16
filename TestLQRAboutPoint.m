q0 = [0 ; 0];
qstar = [pi-1;0];
Q = diag([1 1]);
R = 1;
tf = 5;
umax = 1000;
system = Pendulum(qstar,Q,R);

[t, x] = SimulateLQRSinglePoint(system,q0,qstar,umax,tf );
system.plot(t, x);