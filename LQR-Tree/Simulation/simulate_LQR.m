function [t,x] = simulate_LQR(system,x0,xstar,umax,tf)
    opts = odeset('MaxStep', 0.1, 'RelTol', 1e-4, 'AbsTol', 1e-4);
    [t, x] = ode45(@(t,x) dynamics(system, x,xstar,umax), [0 tf], x0, opts); 
end

function xdot = dynamics(system,x,xstar,umax)
  
    f = system.dynamics;
    u = -inv(system.R)*system.B'*system.S*(x-xstar);
    u(u>umax) = umax;
    u(u<-umax) = -umax;
    xdot = f(x,u);
end

