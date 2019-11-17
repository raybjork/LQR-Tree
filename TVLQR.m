function [K, u] = TVLQR(Q, R, Qf, tf, x_d, u_d, u_max, system)
    L0square = chol(Qf)';

    % ode45() must take L as a vector, so we reshape it
    L0 = reshape(L0square, [2*length(L0square) 1]);
    opts = odeset('MaxStep', 1);
    
    % L must be integrated backwards; we integrate L(tf - t) from 0 to tf
    [t, L] = ode45(@(t,L) dLdtminus(t, L, Q, R, x_d, u_d, system), [0 tf], L0, opts);
    
    % now we flip backwards time to be forwards...
    t = tf - t;
    
    % and rearrange (t, L) to be increasing in time 
    t = flipud(t);
    L = spline(t, flipud(L)');

    K = @(t) gain(system, t, ppval(L, t), R, x_d);
    u = @(t,x) feedback(system, t, x, ppval(L, t), R, x_d, u_d,u_max);
end

function Ldotminus = dLdtminus(t, L, Q, R, x_d, u_d, system)
    % reshape L to a square
    L = reshape(L,[length(L)/2 length(L)/2]);
    
    % compute \dot L
    x = ppval(x_d, t);

    sys_i = system.new(x);

    A = sys_i.A;
    B = sys_i.B;

    dLdt = -.5*(Q*inv(L'))-(A'*L)+.5*(L*L'*B*inv(R)*B'*L);
    
    % set derivative to -\dot L reshaped as a vector
    Ldotminus = -reshape(dLdt,[2*length(dLdt) 1]);
end

function u = feedback(system, t, x, L, R, x_d, u_d,u_max)
    K = gain(system, t, L, R, x_d);
    u =-K*(x-ppval(x_d,t)) + ppval(u_d,t);
    u(u>u_max) = u_max;
    u(u<-u_max)= -u_max;
end

function K = gain(system, t, L, R, x_d)
    sys_i = system.new(ppval(x_d, t));
    B = sys_i.B;
    Lsquare = reshape(L, [length(L)/2 length(L)/2]);
    K = inv(R)*B'*(Lsquare*Lsquare');
end