classdef Pendulum
    properties
        A
        B
        Q
        R
        constants
    end
    
    methods (Access = public)
        function self = Pendulum(A, B)
            self.constants = constants();
        end
        
        function [K, S] = lqr(self, Q, R)
            [K, S] = lqr(self.A ,self.B, Q, R);
        end
    end
       
    methods (Static)
        function dx = f(x, u)
            x(1) = mod(x(1),2*pi);
            c = constants();
            g = c.g;
            m = c.m;
            b = c.b;
            L = c.L;
            dx = [x(2); (-b*x(2) - m*g*L*sin(x(1)) + u)/(m*L^2)];
        end
    end
end

function c = constants()
    c.g = 9.81;
    c.m = 1;
    c.L = 1;
    c.b = 0;
end