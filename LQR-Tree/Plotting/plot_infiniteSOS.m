clear
clc
close all

load('infinite_SOS_data.mat')

sos_color = 1/255*[131 223 255];

N = 300;
[X1,X2] = meshgrid(linspace(0,6,N), linspace(-6,6,N));


hold on
VPLOT = reshape(dmsubs(-V,x,[X1(:) X2(:)]'),size(X1));
[~,h] = contourf(X1,X2,VPLOT,-double(rho)*[1 1]);
colormap(sos_color)
hold on
success = [0,0];
for i = -6:.5:6
    for j = 0:.5:6
        [t,x] = simulate_LQR(system, [j;i],system.qstar,system.u_max,10);
        %plot(x(:,1),x(:,2),'Color','b');
        if(norm(x(end,:)' - system.qstar) < .1)
            if(success(1,1) == 0)
                success = [j,i];
            else
                success = [success; j,i];
            end
        end
        
    end
end

k = boundary(success(:,1),success(:,2));

sim_bound = plot(success(k,1),success(k,2), 'Color' , 'k', 'LineWidth',3);
set(h,'Color',sos_color,'LineWidth',3)

sets = [5,-5;
        2,4.3;
        3,-2;
        3,4];
for i = 1:4
    

    [t,x] = simulate_LQR(system, sets(i,:)',system.qstar,system.u_max,10);
    traj = plot(x(:,1),x(:,2),'Color','b','LineWidth',2);
    start_pt = plot(sets(i,1),sets(i,2),'go','LineWidth',6);
    end_pt = plot(x(end,1),x(end,2),'ro','LineWidth',6);
end

legend([h,sim_bound,traj,start_pt,end_pt],'SOS ROA','Simulated ROA','Sample Trajectories','Start','End')


xlabel('$\theta$','interpreter','latex')
ylabel('$\dot{\theta}$','interpreter','latex')
title('SOS and Simulated ROA for Infinite-Time LQR with Torque Limits');
set(gca,'FontSize',12);

saveas(gcf,'InfiniteTimeSOSFig.png');