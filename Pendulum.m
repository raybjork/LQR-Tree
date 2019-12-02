%% Pendulum
%   class to abstract away dynamics of system in an encapsulated object
classdef Pendulum
    properties
        qstar
        A
        B
        Q
        R
        S
        u_max
        constants
    end

    methods (Access = public)
        %% Constructor
        %   initialize class with a point to linearize about and cost
        %   matrices
        %   can call with three arguments (q, Q, R) or 1 (q)
        function self = Pendulum(varargin)
            if nargin > 0
                self.qstar = varargin{1};
                self.constants = constants();
                [A,B] = linearize(self.qstar);
                self.A = A;
                self.B = B;
                if nargin > 1
                    self.Q = varargin{2};
                    self.R = varargin{3};
                    self.u_max = varargin{4};
                    [~, self.S, ~] = lqr(self.A, self.B, self.Q, self.R);
                end
            end
        end
    end
       
    methods (Static) % Static Methods can be called by [system].[method]
        function system = new(q, u)
            system = Pendulum(q);
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

            % define pendulum object
            L = 1;
            pend = .2;
            theta = x(:,1)';
            px = L*sin(theta);
            py = - L*cos(theta);

            % set bounds of plot
            B = .5;
            range = [-L-B, L+B];

            % initialize video with parameters
            stale = .01;
            tic
            i = 1;

            while i<=numel(t)
                %reset clock
                start = toc;

                % update system
                hold off;
                plot([0, L*sin(theta(i))],[0, -L*cos(theta(i))], ...
                    'k', 'LineWidth',3);
                hold on;
                rectangle('Position',[L * sin(theta(i)) - pend/2, ...
                    -L*cos(theta(i))-pend/2,pend,pend],...
                    'Curvature',[1,1], 'FaceColor','r','EdgeColor', ...
                    'k','LineWidth',3);
                plot(px(1:i), py(1:i), 'g','LineWidth',3);
                axis equal;
                xlim(range);
                ylim(range);
                xlabel('y');
                ylabel('z');
                title(sprintf('Pendulum Trajectory, $t =  %.2f $',t(i)));

                % handle timing of next frame
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
    %% run this section only to get expression for A
    %syms u
    %syms x [2 1]
    %A = jacobian(f(x,u),x); 
    %A = double(subs(A, x, q));
    %B = diff(f(x,u),u);
    %B = double(subs(B, x, q));
    
    %% explicitly derived derivatives
    A = [0,  1; -(981*cos(q(1)))/100, -1];
    B = [0; 1];
end

%% f
%   function to evaluate system dynamics at given state and input
function dx = f(x, u)
    %x(1) = mod(x(1),2*pi); % want to keep theta between 0 and 2pi
    
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
    c.m = 1;
    c.L = 1;
    c.b = 1;
end