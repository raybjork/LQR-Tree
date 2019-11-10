classdef Pendulum
    properties
        A
        B
        Q
        R
        K
        S
        constants
    end
    
    methods (Access = public)
        function self = Pendulum(q,Q,R)
            self.constants = constants();
            [A,B] = linearize(q);
            self.A = A;
            self.B = B;
            [K,S] = lqr(A,B,Q,R);
            self.K = K;
            self.S = S;
            self.Q = Q;
            self.R = R;
            
        end
       
    end
       
    methods (Static)

        
        function dx = Poly_f(x, u)
            %x(1) = mod(x(1),2*pi);
            c = constants();
            g = c.g;
            m = c.m;
            b = c.b;
            L = c.L;
            dx = [x(2); (-b*x(2) - m*g*L*(x(1)+ x(1)^3/6 +x(1)^5/120) + u)/(m*L^2)];
            
        end
    end
end

function [A, B] = linearize(q)
    syms u x1 x2
    x = [x1;x2];
    A = [diff(True_f(x,u),x1),diff(True_f(x,u),x2)];
    A = double(subs(A, x1, q(1)));
    A = double(subs(A, x2, q(2)));
    B = diff(True_f(x,u),u);
    B = double(subs(B, x1, q(1)));
    B = double(subs(B, x2, q(2)));
end

function dx = True_f(x, u)
    %x(1) = mod(x(1),2*pi);
    c = constants();
    g = c.g;
    m = c.m;
    b = c.b;
    L = c.L;
    dx = [x(2); (-b*x(2) - m*g*L*sin(x(1)) + u)/(m*L^2)];
end

function c = constants()
    c.g = 9.81;
    c.m = 1;
    c.L = 1;
    c.b = 1;
end