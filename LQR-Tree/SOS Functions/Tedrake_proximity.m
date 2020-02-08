function [j] = Tedrake_proximity(x_0, x_d,R, system,t_max)
    
    sys = system.new(x_0);
    A = sys.A;
    B = sys.B;
    f = sys.dynamics();
    c = f(x_0,0);
    
    S = get_cost(A,B,c,R,t_max);
    
    j = Inf;
    for t = 0:.001:t_max
        d = get_d(A,c,x_d-x_0,t);
        
        j = min((t+.5*d'*S(t)*d),j);
        
    end
    
    
end


function [S] = get_cost(A,B,c,R,t_max)
    P_0 = A*0;
    P0 = reshape(P_0, [length(P_0)^2 1]);
    opts = odeset('MaxStep', 1);
    
    [t, P] = ode45(@(t,P) Pdot(P, A,B,R), [0 t_max], P0, opts);
    
    
    t = t_max -t;
    t = flipud(t);
    P = spline(t, flipud(P)');

    S = @(t) s(ppval(P,t));
    
    
end

function [P_dot] = Pdot(P,A,B,R)
    P = reshape(P,[sqrt(length(P)) sqrt(length(P))]);
    P_dot = A*P + P*A'+B*inv(R)*B';
    
    P_dot = reshape(P_dot,[length(P_dot)^2 1]);
    
end

function S = s(P)
    Psquare = reshape(P, [sqrt(length(P)) sqrt(length(P))]);
    S = inv(Psquare);
end

function [d] = get_d(A,c,x_d,t_max)
    d = exp(A*t_max)*x_d;%+integral(@(t) exp(A*(t_max - t))*c,0,t_max,'ArrayValued',true);
    d_int = d*0;
    for i = 1:t_max/.01 + 1
      t = (i-1)*.01;
      d_int(:,i) = exp(A*(t_max - t)*c);
    end
    int = .01*cumtrapz(d_int')';
    d = d+int(:,end);
        
end

