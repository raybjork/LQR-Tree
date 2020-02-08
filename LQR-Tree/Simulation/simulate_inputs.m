%% simulate
%   simulate dynamics of system, given controller initial state
function [t, x] = simulate_inputs(system, u, dt, x_0)
    opts = odeset('MaxStep', 0.1, 'RelTol', 1e-4, 'AbsTol', 1e-4);
    [t, x] = ode45(@(t,x) simulation_dynamics(system, t, x, u, dt), ...
        [0 dt * length(u)], x_0, opts); 
end

function dx = simulation_dynamics(system, t, x, u, dt)
    % if state is close to swing-up equilibrium, let LQR take over
    % use cost-to-go as a metric of closeness
    x_star = [0; pi; 0; 0];
    if (x-x_star)'*system.S*(x-x_star) < 1
        u_x = - system.K*(x-x_star);
    else
        % find the right spline i to sample from
        i = ceil(t/dt);
        if i == 0
            i = 1;
        end
        if i >= numel(u)
            i = numel(u) - 1;
        end
        lambda = mod(t,dt)/dt;
        
        % time is between u(i) and u(i+1); select input as r_i(t)
        u_x = u(i)*(1-lambda) + u(i+1)*lambda;
    end
    f = system.dynamics();
    dx = f(x, u_x);
end