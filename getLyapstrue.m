function getLyapstrue(x, system, AB, S,R, Q, dt, u_d, x_d, N)
   V = zeros(1,N)*x(1);
   V_dot = zeros(1,N)*x(1);
   x
   for i = 1:N
       t= (i-1)*dt;
       q_0 = x_d(t);
       [A,B] = AB(t);
       K = R\B'*S(t);
       xhat = x-q_0;
       f = system.dynamics();
       f = system.poly_f(x, u_d(t) - K*xhat);
       V(i) = .5*xhat'*S(t)*xhat;
       S_dot = -(S(t)*A+ A'*S(t) + Q - S(t)*B*inv(R)*B'*S(t));
       V_dot(i) = jacobian(V(i), x) * f + .5 * xhat' * S_dot*xhat;
       eig(double(subs(jacobian(jacobian(V_dot(i), x)', x), x, q_0)))
   end
    
    
    
    
    
    
end