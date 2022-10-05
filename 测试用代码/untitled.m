% clc
% clear
% close all
% h=0.001;
% t=0:h:10;
% Fre=1;
% omega=2*pi*Fre;
% y=abs(sin(omega.*t));
% g=sin(omega.*t);
% figure
% plot(t,y)

% ydot=diff(y)/h;
% ydot(end+1)=ydot(end);

% figure
% plot(t,ydot,'b')

% z=omega.*cos(omega.*t).*(g)./(abs(g));
% % z=omega.*cos(omega.*t);
% hold on
% plot(t,z,'r')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: Shengyi Xu xushengyichn@outlook.com
%Date: 2022-06-16 11:37:30
%LastEditors: Shengyi Xu xushengyichn@outlook.com
%LastEditTime: 2022-06-20 16:46:08
%FilePath: \twindeck_ID\twindeck_ID.m
%Description: ?????????????????


%Copyright (c) 2022 by Shengyi Xu xushengyichn@outlook.com, All Rights Reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

%NTNU笔记本路径
addpath(genpath("C:\Users\shengyix\OneDrive\NAS云同步\Drive\0研究生\有用的代码\HHT-Tutorial-master"))%经验包络法求瞬时频率和振幅必备
addpath("C:\Users\shengyix\OneDrive\NAS云同步\Drive\0研究生\有用的代码\Matlab_PlottingTemplates")%绘图工具
%台式机路径
addpath(genpath("C:\Users\xushe\OneDrive\NAS云同步\Drive\0研究生\有用的代码\HHT-Tutorial-master"))
addpath("C:\Users\xushe\OneDrive\NAS云同步\Drive\0研究生\有用的代码\Matlab_PlottingTemplates")
%% 响应重构
% load a.mat
fs=2000;slength=10;t=1/fs:1/fs:slength;t=t';
told=t;
dt=1/fs;
front_end=0.2;end_end=4;
load xvmd.mat
UP=xvmd;

% % 结构参数
% % Structural parameters
% D=0.667;% deck depth
% m=80;% mass of the segment model
% F0=5.263;% Frequency without wind
% omega0=2*pi*F0;% Circular frequency without wind
% rho=1.225;% density of the air
% U = 6.21; % wind speed
% Zeta0 = 0; % damping ratio without wind
% Fre=5.19865;
% Mass=m;
% h=0.01;
% t=0:h:500;
% t=t';
% u0=0.0001;
% udot0=0;
% P=zeros(1,size(t,1));
% out = polynomial_NB(Fre, Mass, Zeta0, rho, D, U, a, t, P,  u0, udot0);
% % close all
% figure
% plot(out(:,1),out(:,2))
% % hold on 
% % plot(t,UP)
% % legend("measured","calculated")



% %% 提取工况设置


%     UP(:, 1) = out(:,2); %m
% %     DOWN(:, 1) = (D(:, 5) - mean(D(:, 5)) + D(:, 6) - mean(D(:, 6)) + D(:, 7) - mean(D(:, 7)) + D(:, 8) - mean(D(:, 8))) / 4/1000; %

%     AM(:, 2) = UP(:, 1);
% %     AM(:, 3) = DOWN(:, 1);

%     t=out(:, 1);
% fs=1/h;
% [psd_avg, f, psd_plot] = fft_transfer(fs,UP);
% figure
% plot(f, psd_plot)
% % 
% % Wc1=2*1/fs;                                          %下截止频率 1Hz
% % Wc2=2*15/fs;                                          %上截止频率 6Hz
% % [bb,aa]=butter(2,[Wc1, Wc2],'bandpass');  % 二阶的巴特沃斯带通滤波
% % x1=filter(bb,aa,UP);

% % figure
% % plot(t,x1)
% %% 绘制响应曲线
% % The standard values for colors saved in PLOT_STANDARDS() will be accessed from the variable PS
% PS = PLOT_STANDARDS();
% figure(1)
% fig1_comps.fig = gcf;
% hold on
% fig1_comps.p1 = plot(AM(:, 1), AM(:, 2));
% % fig1_comps.p2 = plot(AM(:, 1), AM(:, 3));
% hold off

% %========================================================
% % ADD LABELS, TITLE, LEGEND
% title('Response of the twin girders');
% xlabel('Time/s');
% ylabel('Displacement/m');

% legend([fig1_comps.p1], 'windward');
% legendX = .82; legendY = .87; legendWidth = 0.01; legendHeight = 0.01;
% fig1_comps.legendPosition = [legendX, legendY, legendWidth, legendHeight];
% % If you want the tightest box set width and height values very low matlab automatically sets the tightest box

% %========================================================
% % SET PLOT PROPERTIES
% % Choices for COLORS can be found in ColorPalette.png
% set(fig1_comps.p1, 'LineStyle', '-', 'LineWidth', 2, 'Color', PS.Blue4);
% % set(fig1_comps.p2, 'LineStyle', '-', 'LineWidth', 2, 'Color', PS.MyRed);


% %========================================================
% % INSTANTLY IMPROVE AESTHETICS-most important step
% STANDARDIZE_FIGURE(fig1_comps);


%% 识别气动力参数
%% Identify aerodynamic parameters

% 结构参数
% Structural parameters
% D=0.667;% deck depth
% m=80;% mass of the segment model
% F0=5.263;% Frequency without wind
% omega0=2*pi*F0;% Circular frequency without wind
% rho=1.225;% density of the air
% U = 6.21; % wind speed
% Zeta0 = 0; % damping ratio without wind

D=0.075;% deck depth
m=7.19/0.986;% mass of the segment model
F0=7.88;% Frequency without wind
omega0=2*pi*F0;% Circular frequency without wind
rho=1.22;% density of the air
U =4.81; % wind speed
Zeta0 = 2.34/100; % damping ratio without wind
% Zeta0 = 0; % damping ratio without wind

front_end=0.2;end_end=4;
[ex, frex] = ee(UP, 1 / fs); %经验包络法求瞬时频率和瞬时振幅 % empirical envelope method to find instantaneous frequency and instantaneous amplitude
ex=ex/D;
% front_end = 0; end_end = 0;
slength=t(end);
frex = frex(fs * front_end + 1:fs * (slength - end_end)); %频率边界效应
ex = ex(fs * front_end + 1:fs * (slength - end_end)); %振幅边界效应
UP = UP(fs * front_end + 1:fs * (slength - end_end));
t = t(1:(slength - front_end - end_end) * fs);


% bex=polyfit(t,ex,7);%对振幅曲线进行拟合 TODO：改为 gompertz model 试试
% ex=polyval(bex,t);

% figure
% plot(t, UP)
% hold on
% plot(t, ex, 'r')

[psd_avg, f, psd_plot] = fft_transfer(fs,UP);


omgx = frex * 2 * pi; %瞬时圆频率
bomgx=polyfit(ex,omgx,4);
omgxeq=polyval(bomgx,ex);%瞬时圆频率多项式拟合
% figure;plot(ex,omgx,'g');hold on;plot(ex,omgxeq,'r');title('瞬时频率结果 计算结果(绿)+多项式拟合结果(红)+真实值(蓝)')

epsx=zeros(1,length(ex));epsx=epsx';%瞬时阻尼比
for i=1:length(ex)-1
epsx(i)=log(ex(i)/ex(i+1))./omgx(i)*fs;
end
epsx(length(ex))=epsx(length(ex)-1);
bepsx=polyfit(ex,epsx,4);
epsxeq=polyval(bepsx,ex);%瞬时阻尼比多项式拟合
figure;plot(ex*D,epsx,'g');hold on;plot(ex*D,epsxeq,'r');title('瞬时阻尼结果 计算结果(绿)+多项式拟合结果(红)')

amp=0.0001:0.0001:1;amp=amp';epsxeq=polyval(bepsx,amp);
figure
plot(amp,epsxeq)

for k1=1:length(bepsx)
    c(k1)=bepsx(end-k1+1);%拟合的多项式系数进行倒序排列
end

% a1=-2*(c(1)-Zeta0)*omega0*m/(rho*U*D);
% a2=-3*c(2)*pi*omega0*m/(2*rho*U*D);
% a3=-8*c(3)*omega0*m/(rho*U*D);
% a4=-15*c(4)*pi*omega0*m/(4*rho*U*D);
% a5=-16*c(5)*omega0*m/(rho*U*D);

a1=-2*(c(1)-Zeta0)*omega0*m/(rho*U*D);
a2=-3*c(2)*pi*omega0*m/(2*rho*U*D);
a3=-8*c(3)*omega0*m/(rho*U*D);
a4=-15*c(4)*pi*omega0*m/(4*rho*U*D);
a5=-16*c(5)*omega0*m/(rho*U*D);
tt=told;
% told=0:(1/1000000):10;
% told=told';
a=[a1 a2 a3 a4 a5];
% a=[a1 a2*0.99 a3*0.98 a4 a5];
P=zeros(1,size(told,1));

amp=0.0001:0.0001:1;amp=amp';
epsall=-rho.*U.*D.*(a(1)+4.*a(2).*amp./3./pi+a(3).*amp.^2/4+8.*a(4).*amp.^3/15/pi+a(5).*amp.^4/8)./2./omega0./m+Zeta0;
figure
plot(amp*D,epsall)
hold on
plot(ex*D,epsx,'g')

u0=xvmd(1000);
udot0=(xvmd(1001)-xvmd(1000))/dt;
tt=told(1000:end)-told(1000);
dis=xvmd(1000:end);

% u0=xvmd(1);
% udot0=(xvmd(2)-xvmd(1))/dt;

% u0=0.001;
% udot0=0;

Fre=F0;
Mass=m;
out = polynomial_NB(Fre, Mass, Zeta0, rho, D, U, a, told, P,  u0, udot0);
% close all
figure


plot(out(:,1),out(:,2),'b')
hold on
plot(tt,dis,'r')
xlabel("时间/s")
ylabel("振幅/m")
