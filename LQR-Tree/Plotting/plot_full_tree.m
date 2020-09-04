clear
clc
close all

sos_color = 1/255*[131 223 255];

% load data from a saved *.mat file
%load examplefile.mat

for i = 1:N
    hold on
    plot_size = 10;
    density = 500;
    [X1,X2] = meshgrid(linspace(-plot_size,plot_size,density), linspace(-plot_size,plot_size,density));
    VPLOT = reshape(dmsubs(-V(i),x,[X1(:) X2(:)]'),size(X1));
    [~,sos] = contourf(X1,X2,VPLOT,-double(rho(i))*[1 1]);
    colormap(sos_color)

    set(sos,'Color',sos_color,'LineWidth',3)
end

% load data from a saved *.mat file
%load examplefile.mat

for i = 1:N
    
    plot_size = 10;
    density = 500;
    [X1,X2] = meshgrid(linspace(-plot_size,plot_size,density), linspace(-plot_size,plot_size,density));
    VPLOT = reshape(dmsubs(-V(i),x,[X1(:) X2(:)]'),size(X1));
    [~,sos] = contourf(X1,X2,VPLOT,-double(rho(i))*[1 1]);
    colormap(sos_color)

    set(sos,'Color',sos_color,'LineWidth',3)
end

% load data from a saved *.mat file
%load examplefile.mat

for i = 1:N
    
    plot_size = 10;
    density = 500;
    [X1,X2] = meshgrid(linspace(-plot_size,plot_size,density), linspace(-plot_size,plot_size,density));
    VPLOT = reshape(dmsubs(-V(i),x,[X1(:) X2(:)]'),size(X1));
    [~,sos] = contourf(X1,X2,VPLOT,-double(rho(i))*[1 1]);
    colormap(sos_color)

    set(sos,'Color',sos_color,'LineWidth',3)
end
% load data from a saved *.mat file
%load examplefile.mat
state = x_d(0:dt/5:dt*(N-1));
b = plot(state(1,:),state(2,:),'Color','b','LineWidth',2);

% load data from a saved *.mat file
%load examplefile.mat
state = x_d(0:dt/5:dt*(N-1));
g = plot(state(1,:),state(2,:),'Color','g','LineWidth',2);

% load data from a saved *.mat file
%load examplefile.mat
state = x_d(0:dt/5:dt*(N-1));
k = plot(state(1,:),state(2,:),'Color','k','LineWidth',2);

ylim([-6,6])
xlim([-3,9])

legend([sos,b,g,k],'Trajectory ROA','Trajectory 1', 'Trajectroy 2', 'Trajectory 3');


title('LQR-Tree example with 3 trajectories');
xlabel('$\theta$','interpreter','latex')
ylabel('$\dot{\theta}$','interpreter','latex')
set(gca,'FontSize',12);

saveas(gcf,'LQRTree3Traj.png');