clc
clear
close all
load('results.mat')


% scatter_handle = findobj(fig_handle, 'Type', 'scatter3');
% xdata = get(scatter_handle, 'XData');
% ydata = get(scatter_handle, 'YData');

xdata=cases(:,1);


ydata=cases(:,2);


zdata=results;


% scatter3(x,y,z)

% xdata=x;
% ydata=y;
% zdata=z;

scatter3(xdata, ydata, zdata, [], zdata, 'filled');

xlim([min(xdata), max(xdata)]); % replace xmin and xmax with the desired limits for the x axis
ylim([min(ydata), max(ydata)]); % replace ymin and ymax with the desired limits for the y axis
zlim([min(zdata), max(zdata)]); % replace zmin and zmax with the desired limits for the z axis

xlabel("频率（Hz）")
ylabel("位置（m）")
zlabel("振幅（m）")

% [X, Y] = meshgrid(linspace(min(xdata), max(xdata), 500), linspace(min(ydata), max(ydata), 500));
% Z = griddata(xdata, ydata, zdata, X, Y, 'linear');
% 
% surf(X, Y, Z);