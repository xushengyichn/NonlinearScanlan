clc
clear
close all
load('damping_plot.mat')
FTMD1=0.817503524129623;
FTMD2=FTMD2_all/FTMD1;
figure
surf(U_temps*maxphi1, FTMD2, bridge_damping_grid, 'FaceLighting','gouraud',...
    'MeshStyle','column',...
    'SpecularColorReflectance',0,...
    'SpecularExponent',5,...
    'SpecularStrength',0.2,...
    'DiffuseStrength',1,...
    'AmbientStrength',0.4,...
    'AlignVertexCenters','on',...
    'LineWidth',0.2,...
    'FaceAlpha',1,...
    'FaceColor','r',...
    'EdgeAlpha',0.2);
hold on
shading interp % 美化1
camlight       % 美化2
lighting phong % 美化3
surf(U_temps2_single*maxphi1, FTMD2, bridge_damping_grid_single,...
    'AmbientStrength',0.4,...
    'FaceColor',[1 1 1],...
    'AlignVertexCenters','on',...
    'FaceAlpha',0.6,...
    'LineWidth',0.2,...
    'EdgeAlpha',1);
shading interp % 美化1

colormap(winter)

% 
% % plot parameter -------------------------------
% xmin= -10;xmax = 15;
% ymin= -10;ymax = 15;
% zmin=  0;zmax = 15;
% 
% % view point ----------------------------------
% vx = -50; vy = 20;
% 
% 
% fig_pos = [0.09 0.15 0.53 0.75];
% axis_pos= [xmin xmax ymin ymax zmin zmax];
% 
% 
% 
% 
% % plot 3D points ------------------------------
% h=figure(1);
% set(h, 'Position', [100, 100, 800, 500]);
% 
% % plot zeros plane ------------------------------
% ax1 = axes;
% set(ax1,'position',fig_pos );
% h_1 = surf(U_temps*maxphi1, FTMD2, bridge_damping_grid, 'FaceLighting','gouraud',...
%     'MeshStyle','column',...
%     'SpecularColorReflectance',0,...
%     'SpecularExponent',5,...
%     'SpecularStrength',0.2,...
%     'DiffuseStrength',1,...
%     'AmbientStrength',0.4,...
%     'AlignVertexCenters','on',...
%     'LineWidth',0.2,...
%     'FaceAlpha',1,...
%     'FaceColor','r',...
%     'EdgeAlpha',0.2);
% box on
% hidden off
% % axis(axis_pos );
% view([vx,vy])
% colormap(ax1,summer);
% xlabel('IamX')
% ylabel('IamY')
% zlabel('IamZ')
% title('IamTitle')
% 
% 
% 
% % % plot x plane ------------------------------
% % ax2 = axes;
% % set(ax2,'position',fig_pos);
% % h_2 = surf(U_temps2_single*maxphi1, FTMD2, bridge_damping_grid_single,...
% %     'AmbientStrength',0.4,...
% %     'FaceColor',[1 1 1],...
% %     'AlignVertexCenters','on',...
% %     'FaceAlpha',0.6,...
% %     'LineWidth',0.2,...
% %     'EdgeAlpha',1);
% % box off
% % axis off
% % hidden off
% % % axis(axis_pos );
% % view([vx,vy])
% % colormap(ax2,'pink');
% % 
% % % set colorbar --------------------------------
% % cb1 = colorbar(ax1,'Position',[.65  .15 .05 .6]);
% % cb2 = colorbar(ax2,'Position',[.74  .15 .05 .6]);
% % % cb3 = colorbar(ax3,'Position',[.83   .15 .05 .6]);
% % % cb3 = colorbar(ax4,'Position',[.92   .15 .05 .6]);
% % 
% % % set fontsize --------------------------------
% % set(ax1,'fontsize',mm_fz)
% % set(ax2,'fontsize',mm_fz)
% % % set(ax3,'fontsize',mm_fz)
% % % set(ax4,'fontsize',mm_fz)


a=U_temps*maxphi1
b= FTMD2
c=bridge_damping_grid
bridge_damping_grid_single=griddata(variables2(:,1),variables2(:,2),singleplotdata_min_rep,U_temps2_single,FTMD2_all2);
figure
variables=[A(:) B(:) C(:)]