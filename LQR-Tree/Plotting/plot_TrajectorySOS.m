clear
clc
close all

% load data from a saved *.mat file
%load examplefile.mat
N = 300;
[X1,X2] = meshgrid(linspace(0,6,N), linspace(-6,6,N));

VPLOT = reshape(dmsubs(-V,x,[X1(:) X2(:)]'),size(X1));
[~,h] = contour(X1,X2,VPLOT,-double(rho)*[1 1]);
set(h,'Color','k','LineWidth',3)

 load trajectorySOSdata.mat
% 
 sos_color = 1/255*[131 223 255];
 hold on

 title('SOS ROA For Nominal Trajectory: Simultaneous Generation','interpreter','latex');
 xlabel('$\theta$','interpreter','latex')
 ylabel('$\dot{\theta}$','interpreter','latex')
 set(gca,'FontSize',12);

for i = 1:N
    
    plot_size = 6;
    density = 300;
    [X1,X2] = meshgrid(linspace(-plot_size,plot_size,density), linspace(-plot_size,plot_size,density));
    VPLOT = reshape(dmsubs(-V(i),x,[X1(:) X2(:)]'),size(X1));
    [~,sos] = contourf(X1,X2,VPLOT,-double(rho(i))*[1 1]);
    colormap(sos_color)

    set(sos,'Color',sos_color*.9,'LineWidth',3)
    
    
    traj = plot(state(1,:),state(2,:),'Color','b','LineWidth',2);
    hold on
    frame(i) = getframe(gcf);
    
    l = legend([h,sos,traj],'Infinite-Time LQR ROA','SOS ROA','Nominal Trajectory');
    set(l,'Location','southwest');
    
end


video = VideoWriter('TrajectorySOS.avi', 'Uncompressed AVI');
video.FrameRate = 5;
open(video)
writeVideo(video, frame);
close(video);



