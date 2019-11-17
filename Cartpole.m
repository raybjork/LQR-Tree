        methods (Access = public)
        function self = Cartpole(A, B)
            self.constants = constants();

            self.A = [0 0 1 0;
                  0 0 0 1;
                  0 self.constants.g*self.constants.mp/self.constants.mc 0 0;
                  0 self.constants.g*(self.constants.mc+self.constants.mp)/(self.constants.L*self.constants.mc) 0 0];
            self.B = [0;0; 1/self.constants.mc; 1/(self.constants.L*self.constants.mc)];

        end
        
        function [K, S] = lqr(self, Q, R)
            [K, S] = lqr(self.A ,self.B, Q, R);
        end
    end
       
        
        function dx = f( x, K)
            x(2) = mod(x(2),2*pi);
            c = constants();
            g = c.g;
            mc = c.mc;
            mp = c.mp;
            L = c.L;
            M = [mc + mp, mp*L*cos(x(2));
                 mp*L*cos(x(2)), mp*L^2];
            C = [-mp*L*sin(x(2))*x(4)^2;
                 mp*g*L*sin(x(2))];
            B = [1;
                 0];

            u = - K*(x-[0;pi;0;0]);
            dx = [x(3:4);
                  M \ (B*u - C)];
        end
    end

%% Cartpole
%   class to abstract away dynamics of system in an encapsulated object
classdef Cartpole
    properties
        qstar
        A
        B
        Q
        R
        K
        S
        constants
    end

    methods (Access = public)
        %% Constructor
        %   initialize class with a point to linearize about and cost
        %   matrices
        %   can call with three arguments (q, Q, R) or 1 (q)
        function self = Cartpole(varargin)
            self.qstar = varargin{1};
            self.constants = constants();
            [A,B] = linearize(self.qstar);
            self.A = A;
            self.B = B;
            if (nargin == 3)
                self.Q = varargin{2};
                self.R = varargin{3};  
                [K,S] = lqr(A, B, self.Q, self.R);
                self.K = K;
                self.S = S;
            end
        end
    end
       
    methods (Static) % Static Methods can be called by [system].[method]
        function system = new(q)
            system = Cartpole(q);
        end
        
        %% dynamics
        %   return a handle to the function encapsulating the dynamics
        function handle = dynamics()
            handle = @f;
        end
        
        %% poly_f
        %   function to evaluate polynomial approximation of system dynamics
        %   at a given state and input
        function dx = poly_f(x, u)
            %x(1) = mod(x(1),2*pi); want to keep theta between 0 and 2pi

            % unpack constants
            c = constants();
            g = c.g;
            m = c.m;
            b = c.b;
            L = c.L;

            dx = [x(2); (-b*x(2) - m*g*L*(x(1)+x(1)^3/6 + x(1)^5/120) + u)/(m*L^2)];
        end
        
        %% plot
        %   plot trajectory of the system
        function plot(t, x)
            % set text display to Latex
            set(groot,'defaulttextinterpreter','latex');
            set(groot, 'DefaultLegendInterpreter', 'latex')
            
            y_d = 0*t;
            z_d = 0*t;

            th = x(:,2)';
            x = x(:,1)';

            buffer = 2;
            xrange = [min(x) - buffer, max(x) + buffer];
            yrange = [-1, 3];
            tic

            c = constants();
            L = c.L;
            
            h = .2;
            w = .4;
            pend = .1;
            pennblue = [1,37,110]/256;
            pennred = [149,0,26]/256;
            px = x + L*sin(th);
            py = - L*cos(th);


            stale = .01;
            tic

            i = 1;

            while i<=numel(t)
                start = toc;
                hold off;

                plot(4*xrange,[0 0], 'k', 'LineWidth',3)
                hold on;
                rectangle('Position',[x(i)-w/2, -h/2, w, h],'FaceColor',pennblue,'EdgeColor','k',...
                'LineWidth',3)

                plot([x(i), x(i) + L*sin(th(i))],[0, -L*cos(th(i))], 'k', 'LineWidth',3);

                rectangle('Position',[x(i) + L*sin(th(i))-pend/2,-L*cos(th(i))-pend/2,pend,pend],...
                    'Curvature',[1,1], 'FaceColor',pennred,'EdgeColor','k','LineWidth',3);
                plot(px(1:i), py(1:i), 'g','LineWidth',3);
                axis equal;
                %xlim([x(i) - 2, x(i) + 2]);
                xlim(xrange);
                ylim(yrange);
                %legend('x_d(t)','x(t)');
                xlabel('y');
                ylabel('z');
                titl = sprintf('Cart-Pole Trajectory, $t =  %.2f $',t(i));
                title(titl,'Interpreter','latex');


                compu = toc - start;
                stale_i = max(stale,compu*2);
                next_i = find(t >= start + stale_i);
                if numel(next_i) < 1
                    if i < numel(t)
                        i = numel(t);
                    else
                        break;
                    end
                else
                    i = next_i(1);
                end
                pause(t(i) - toc);

            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Private Methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% linearize
%   linearize system around a given state
function [A, B] = linearize(q)
    syms u
    syms x [2 1]
    %A = jacobian(f(x,u),x); % run this only to get expression for A
    A = [0,  1; -(981*cos(x1))/100, -1];
    A = double(subs(A, x, q));
    %B = diff(f(x,u),u); % run this only to get expression for B
    B = [0; 1];
    B = double(subs(B, x, q));
end

%% f
%   function to evaluate system dynamics at given state and input
function dx = f(x, u)
    %x(1) = mod(x(1),2*pi);% want to keep theta between 0 and 2pi
    
    % unpack constants
    c = constants();
    g = c.g;
    m = c.m;
    b = c.b;
    L = c.L;
    
    dx = [x(2); (-b*x(2) - m*g*L*sin(x(1)) + u)/(m*L^2)];
end

%% constants
%   function to pack constants into struct
function c = constants()
    c.g = 9.81;
    c.mp = 1;
    c.mc = 1;
    c.L = 1;
end