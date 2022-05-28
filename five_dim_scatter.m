% +------------------------------------------------------+
% |                5D Data Visualization                 |
% |              with MATLAB Implementation              | 
% |                                                      |
% | Author: M.Sc. Eng. Hristo Zhivomirov        12/13/14 | 
% +------------------------------------------------------+

clear, clc, close all

% form the axes
x = -1:0.1:1;	% first dimension indipendent variable
y = -1:0.1:1;	% second dimension indipendent variable
z = -1:0.1:1;	% third dimension indipendent variable
t = 0:0.05:1;   % fourth dimension indipendent variable
[X, Y, Z] = meshgrid(x, y, z);

% form the data matrix - it is the fifth dimension
% the data could be imported from a file or could be generated via equation  
load 5D_data

%% organize the visualization
figure(1) 
hScatter = scatter3(X(:), Y(:), Z(:), 'filled');
colormap jet
grid3 on

set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
xlabel('X')
ylabel('Y')
zlabel('Z')

colorbar
caxis([min(Data(:)) max(Data(:))])

% cycle through the fourth dimension independent variable
for k = 1:length(t)	
    
    % visualize the fifth dimension via the marker color
    C = reshape(Data(:, :, :, k), length(x)*length(y)*length(z), 1);  
       
    % update the plot
    set(hScatter, 'CData', C, 'SizeData', 25)  
    title(['Data = \it{f} \rm(X, Y, Z, t) @ t = ' num2str(t(k))])
    drawnow
    
    % pause for a while
    pause(1)   
    
end