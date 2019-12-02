function [S, AB, u] = TVLQR(x_d, u_d, tf, system)
    L0square = chol(system.S)';

    % ode45() must take L as a vector, so we reshape it
    L0 = reshape(L0square, [length(L0square)^2 1]);
    opts = odeset('MaxStep', 1);
    
    % L must be integrated backwards; we integrate L(tf - t) from 0 to tf
    [t, L] = ode45(@(t,L) Ldot(t, L, x_d, u_d, system), [0 tf], L0, opts);
    
    % now we flip backwards time to be forwards...
    t = tf - t;
    
    % and rearrange (t, L) to be increasing in time 
    t = flipud(t);
    L = spline(t, flipud(L)');

    S = @(t) s(ppval(L,t));
    AB = @(t) ab(system, x_d(t), u_d(t));
    u = @(t,x) feedback(t, x, ppval(L, t), x_d, u_d, system);
end

function dLdt = Ldot(t, L, x_d, u_d, system)
    % reshape L to a square
    L = reshape(L,[sqrt(length(L)) sqrt(length(L))]);

    sys_i = system.new(x_d(t), u_d(t));

    A = sys_i.A;
    B = sys_i.B;

    dLdt = -.5*(system.Q*inv(L'))-(A'*L)+.5*(L*L'*B*inv(system.R)*B'*L);
    
    % set derivative to -\dot L reshaped as a vector
    dLdt = -reshape(dLdt,[length(dLdt)^2 1]);
end

function S = s(L)
    Lsquare = reshape(L, [sqrt(length(L)) sqrt(length(L))]);
    S = (Lsquare*Lsquare');
end

function [A, B] = ab(system, q, u)
    sys_i = system.new(q, u);
    A = sys_i.A;
    B = sys_i.B;
end

function u = feedback(t, x, L, x_d, u_d, system)
    if (x - system.qstar)' * system.S * (x - system.qstar) <3 % change me
        u = -inv(system.R) * system.B' * system.S * (x - system.qstar);
    else
        [~, B] = ab(system, x_d(t), u_d(t));
        S = s(L);
        u =-inv(system.R) * B' * S * (x - x_d(t)) + u_d(t);
    end
    u(u > system.u_max) = system.u_max;
    u(u < -system.u_max)= -system.u_max;
end
