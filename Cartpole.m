classdef Cartpole
    properties
        A
        B
        Q
        R
        constants
    end
    
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
       
    methods (Static)
        function [t, x] = simulate(K, tf, x_0)
            opts = odeset('MaxStep', 0.1,'RelTol',1e-4,'AbsTol',1e-4);
            [t, x] = ode45(@(t,x) f(t, x, K), [0 tf], x_0, opts);
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

        function plot(t, x)
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


function c = constants()
    c.g = 9.81;
    c.mp = 1;
    c.mc = 1;
    c.L = 1;
end