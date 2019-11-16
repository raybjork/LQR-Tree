function [t, x] = simulate_TVLQR(system, u, tf, x0, umax)
    opts = odeset('MaxStep', 0.1, 'RelTol', 1e-4, 'AbsTol', 1e-4);
    [t, x] = ode45(@(t,x) dynamics(t, x, system, u, umax), [0 tf], x0, opts); 
end

function xdot = dynamics(t, x, system, u, umax)
    u_x = u(t, x);
    if u_x > umax
        u_x = umax;
    elseif u_x < -umax
        u_x = -umax;
    end
    f = system.dynamics();
    xdot = f(x, u_x);
end