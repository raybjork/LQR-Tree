q0 = [pi-1 ; 0];
qstar = [pi;0];
Q = diag([1 1]);
R = 1;
tf = 5;
umax = 20;
system = Pendulum(qstar,Q,R);

[t, x] = SimulateLQRSinglePoint(system,q0,qstar,umax,tf );
system.plot(t, x);