%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%Date: 2022-11-17 14:05:14
%LastEditors: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%LastEditTime: 2022-11-18 10:51:20
%FilePath: \NonlinearScanlan\20221114一阶模态二TMD结果穷举\eigenvalue_value_of_2tmd_system.m
%Description: 求解2TMD系统的特征值，着重研究负阻尼情况下的TMD最优参数
%
%Copyright (c) 2022 by xushengyichn 54436848+xushengyichn@users.noreply.github.com, All Rights Reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

% 2TMD系统参数
m = 1;
k = 27.449810272952913;
omega = sqrt(k / m);
fs = omega / (2 * pi);
xi_sys = 0.003;
% xi_sys=0.05;
c = 2 * m * omega * xi_sys;

mu1 = 0.02;
mtmd1 = mu1 * m;
f1 = 1 / (1 + mu1) * fs;
xi1 = sqrt(3 * mu1 / 8 / (1 + mu1));
omega1 = 2 * pi * f1;
ctmd1 = 2 * mtmd1 * omega1 * xi1;
ktmd1 = mtmd1 * omega1^2;

% 1TMD系统矩阵
M = [m 0; 0 mtmd1];
C = [c + ctmd1 -ctmd1; -ctmd1 ctmd1];
K = [k + ktmd1 -ktmd1; -ktmd1 ktmd1];
Mode = Complex_Eigenvalue_Analysis(M, C, K);

disp(Mode)
singledamp = min(Mode.("Damping ratio"));

%% 寻找一个TMD的最优参数
xi1s = 0.05:0.01:0.3;
ftmd1s = fs * 0.5:0.01:fs * 2;
[Xi1s, Ftmd1s] = ndgrid(xi1s, ftmd1s);
variables = [Xi1s(:) Ftmd1s(:)];

Damping_ratio = zeros(size(variables, 1), 1);

for k1 = 1:size(variables, 1)
    ctmd1 = 2 * mtmd1 * omega1 * variables(k1, 1);
    omega1 = 2 * pi * variables(k1, 2);
    ktmd1 = mtmd1 * omega1^2;
    M = [m 0; 0 mtmd1];
    C = [c + ctmd1 -ctmd1; -ctmd1 ctmd1];
    K = [k + ktmd1 -ktmd1; -ktmd1 ktmd1];
    Mode = Complex_Eigenvalue_Analysis(M, C, K);
    Damping_ratio(k1, 1) = min(Mode.("Damping ratio"));
end

Damping_ratio_grid = griddata(variables(:, 1), variables(:, 2), Damping_ratio(:, 1), Xi1s, Ftmd1s);
surf(Xi1s, Ftmd1s, Damping_ratio_grid)

maxvalue = max(max(Damping_ratio_grid))
[a, b] = find(Damping_ratio_grid == maxvalue)
f1_opt = Ftmd1s(a, b)
xi1_opt = Xi1s(a, b)

f1 = f1_opt;
xi1 = xi1_opt;
omega1 = 2 * pi * f1;
ctmd1 = 2 * mtmd1 * omega1 * xi1;
ktmd1 = mtmd1 * omega1^2;

% 1TMD系统矩阵
M = [m 0; 0 mtmd1];
C = [c + ctmd1 -ctmd1; -ctmd1 ctmd1];
K = [k + ktmd1 -ktmd1; -ktmd1 ktmd1];
Mode = Complex_Eigenvalue_Analysis(M, C, K);

disp(Mode)

%% 系统频率变化情况下TMD的控制效果
m = 1;
k = 27.449810272952913;

ks = k * 0.7:0.01:k * 1.4;

for k1 = 1:length(ks)
    omega_temp = sqrt(ks(k1) / m);
    fs = omega_temp / (2 * pi);
    c = 2 * m * omega_temp * xi_sys;
    k_temp = ks(k1);
    M = [m 0; 0 mtmd1];
    C = [c + ctmd1 -ctmd1; -ctmd1 ctmd1];
    K = [k_temp + ktmd1 -ktmd1; -ktmd1 ktmd1];
    Mode = Complex_Eigenvalue_Analysis(M, C, K);
    Damping_ratio_change_sys_frequency(k1) = min(Mode.("Damping ratio"));
end

figure
plot(ks / k, Damping_ratio_change_sys_frequency)

%% 非模态坐标

mu2 = 0.02;
mtmd2 = mu2 * m;
f2 = 1 / (1 + mu2) * fs;
xi2 = sqrt(3 * mu2 / 8 / (1 + mu2));
omega2 = 2 * pi * f2;
ctmd2 = 2 * mtmd2 * omega2 * xi2;
ktmd2 = mtmd2 * omega2^2;

xi2s = xi1;
ftmd2s = f1 * 0.7:0.01:f1 * 1.4;
[Xi2s, Ftmd2s] = ndgrid(xi2s, ftmd2s);
variables2 = [Xi2s(:) Ftmd2s(:)];

Damping_ratio2 = zeros(size(variables2, 1), 1);

for k1 = 1:size(variables2, 1)
    ctmd2 = 2 * mtmd2 * omega2 * variables2(k1, 1);
    omega2 = 2 * pi * variables2(k1, 2);
    ktmd2 = mtmd2 * omega2^2;
    M = [m 0 0; 0 mtmd1 0; 0 0 mtmd2];
    C = [c + ctmd1 + ctmd2 -ctmd1 -ctmd2; -ctmd1 ctmd2 0; -ctmd2 0 ctmd2];
    K = [k + ktmd1 + ktmd2 -ktmd1 -ktmd2; -ktmd1 ktmd2 0; -ktmd2 0 ktmd2];
    Mode = Complex_Eigenvalue_Analysis(M, C, K);
    Damping_ratio2(k1, 1) = min(Mode.("Damping ratio"));
    Frequency2(1:3, k1) = Mode.Frequency;

end

single = repmat(singledamp, length(ftmd2s / f1), 1);
figure
plot(ftmd2s / f1, Damping_ratio2)
hold on
plot(ftmd2s / f1, single)

% figure
% Damping_ratio_grid2=griddata(variables2(:,1),variables2(:,2),Damping_ratio2(:,1),Xi2s,Ftmd2s);
% surf(Xi2s/xi1,Ftmd2s/f1,Damping_ratio_grid2)
%
% figure
% hold on
% Frequency2_grid_1=griddata(variables2(:,1),variables2(:,2),Frequency2(1,:),Xi2s,Ftmd2s);
% Frequency2_grid_2=griddata(variables2(:,1),variables2(:,2),Frequency2(2,:),Xi2s,Ftmd2s);
% Frequency2_grid_3=griddata(variables2(:,1),variables2(:,2),Frequency2(3,:),Xi2s,Ftmd2s);
% surf(Xi2s/xi1,Ftmd2s/f1,Frequency2_grid_1)
% surf(Xi2s/xi1,Ftmd2s/f1,Frequency2_grid_2)
% surf(Xi2s/xi1,Ftmd2s/f1,Frequency2_grid_3)
%
%
%% 模态坐标
clc
clear

AD = 0; %气动阻尼-0.8108
% AD=-0.8108;%气动阻尼-0.8108
m = 1;
k = 27.449810272952913;
omega = sqrt(k / m);
fs = omega / (2 * pi);
xi_sys = 0.003;
c = 2 * m * omega * xi_sys;
phi1max = 4.765120284718218e-04;
mu1 = 0.02;
mtmd1 = mu1 / phi1max^2;

f1=1/(1+mu1)*fs;
xi1=sqrt(3*mu1/8/(1+mu1));
% f1 = 0.826927;
% xi1 = 0.07;
omega1 = 2 * pi * f1;
ctmd1 = 2 * mtmd1 * omega1 * xi1;
ktmd1 = mtmd1 * omega1^2;

% 1TMD系统矩阵
M = [m 0; 0 mtmd1];
C = [c + ctmd1 * phi1max^2 + AD -ctmd1 * phi1max; -ctmd1 * phi1max ctmd1];
K = [k + ktmd1 * phi1max^2 -ktmd1 * phi1max; -ktmd1 * phi1max ktmd1];
Mode = Complex_Eigenvalue_Analysis(M, C, K);

disp(Mode)
singledamp = min(Mode.("Damping ratio"));
singlefreq = Mode.Frequency;
disp(singlefreq)

%% 施加第二个TMD

m = 1;
k = 27.449810272952913;
omega = sqrt(k / m);
fs = omega / (2 * pi);

c = 2 * m * omega * xi_sys;
phi1max = 4.765120284718218e-04;
mu1 = 0.02;
mtmd1 = mu1 / phi1max^2;

f1 = 1 / (1 + mu1) * fs;
xi1 = sqrt(3 * mu1 / 8 / (1 + mu1));
omega1 = 2 * pi * f1;
ctmd1 = 2 * mtmd1 * omega1 * xi1;
ktmd1 = mtmd1 * omega1^2;

phi1max = 4.765120284718218e-04;
mu2 = 0.02;
mtmd2 = mu2 / phi1max^2;
mu1 = 0.02;
mtmd1 = mu1 / phi1max^2;

ftmd2s = f1 * 0.7:0.01:f1 * 1.4;
ftmd2s = ftmd2s
xi2 = xi1;
% omega2=2*pi*ftmd2s;
% ctmd2=2*mtmd2*omega2*xi2;
% ktmd2=mtmd2*omega2.^2;

% for k1 = 1:1
for k1 = 1:length(ftmd2s)
    omega2 = 2 * pi * ftmd2s(k1);
    ctmd2 = 2 * mtmd2 * omega2 * xi2;

    ktmd2 = mtmd2 * omega2^2;
    M = [m 0 0; 0 mtmd1 0; 0 0 mtmd2];
    C = [c + ctmd1 * phi1max^2 + ctmd2 * phi1max^2 + AD -ctmd1 * phi1max -ctmd2 * phi1max; -ctmd1 * phi1max ctmd1 0; -ctmd2 * phi1max 0 ctmd2];
    K = [k + ktmd1 * phi1max^2 + ktmd2 * phi1max^2 -ktmd1 * phi1max -ktmd2 * phi1max; -ktmd1 * phi1max ktmd1 0; -ktmd2 * phi1max 0 ktmd2];
    Mode = Complex_Eigenvalue_Analysis(M, C, K);
    Damping_ratio2(k1, 1) = min(Mode.("Damping ratio"));
    Frequency2(1:3, k1) = Mode.Frequency;
end

single = repmat(singledamp, length(ftmd2s / f1), 1);
figure
plot(ftmd2s / f1, Damping_ratio2)
hold on
plot(ftmd2s / f1, single)

figure
plot(ftmd2s / f1, Frequency2(1, :))
hold on
plot(ftmd2s / f1, Frequency2(2, :))
plot(ftmd2s / f1, Frequency2(3, :))

%% 三维绘图
f1 = 1 / (1 + mu1) * fs;
xi1 = sqrt(3 * mu1 / 8 / (1 + mu1));
%     f1=0.826927;%针对阻尼比为-0.06的情况
%     xi1=0.07;
omega1 = 2 * pi * f1;
ctmd1 = 2 * mtmd1 * omega1 * xi1;
ktmd1 = mtmd1 * omega1^2;
AD = 0; %气动阻尼-0.8108
% AD=-0.8108;%气动阻尼-0.8108
xi_syss = -0.1:0.01:0.1;
ftmd2s = f1 * 0.7:0.01:f1 * 1.4;

[Xi_syss, Ftmd2s] = ndgrid(xi_syss, ftmd2s);
variables3 = [Xi_syss(:), Ftmd2s(:)];

for k2 = 1:size(variables3, 1)
    % for k2=1
    m = 1;
    k = 27.449810272952913;
    omega = sqrt(k / m);
    fs = omega / (2 * pi);
    xi_sys = variables3(k2, 1);
    c = 2 * m * omega * xi_sys;
    phi1max = 4.765120284718218e-04;
    mu1 = 0.02;
    mtmd1 = mu1 / phi1max^2;

    % 1TMD系统矩阵
    M = [m 0; 0 mtmd1];
    C = [c + ctmd1 * phi1max^2 + AD -ctmd1 * phi1max; -ctmd1 * phi1max ctmd1];
    K = [k + ktmd1 * phi1max^2 -ktmd1 * phi1max; -ktmd1 * phi1max ktmd1];
    Mode = Complex_Eigenvalue_Analysis(M, C, K);

    %     disp(Mode)
    Damping_ratio1(k2) = min(Mode.("Damping ratio"));

    m = 1;
    k = 27.449810272952913;
    omega = sqrt(k / m);
    fs = omega / (2 * pi);

    c = 2 * m * omega * xi_sys;
    phi1max = 4.765120284718218e-04;
    mu1 = 0.02;
    mtmd1 = mu1 / phi1max^2;

    f1 = 1 / (1 + mu1) * fs;
    xi1 = sqrt(3 * mu1 / 8 / (1 + mu1));
    omega1 = 2 * pi * f1;
    ctmd1 = 2 * mtmd1 * omega1 * xi1;
    ktmd1 = mtmd1 * omega1^2;

    phi1max = 4.765120284718218e-04;
    mu2 = 0.02;
    mtmd2 = mu2 / phi1max^2;
    mu1 = 0.02;
    mtmd1 = mu1 / phi1max^2;

    ftmd2s = variables3(k2, 2);
    xi2 = xi1;

    omega2 = 2 * pi * ftmd2s;
    ctmd2 = 2 * mtmd2 * omega2 * xi2;

    ktmd2 = mtmd2 * omega2^2;
    M = [m 0 0; 0 mtmd1 0; 0 0 mtmd2];
    C = [c + ctmd1 * phi1max^2 + ctmd2 * phi1max^2 + AD -ctmd1 * phi1max -ctmd2 * phi1max; -ctmd1 * phi1max ctmd1 0; -ctmd2 * phi1max 0 ctmd2];
    K = [k + ktmd1 * phi1max^2 + ktmd2 * phi1max^2 -ktmd1 * phi1max -ktmd2 * phi1max; -ktmd1 * phi1max ktmd1 0; -ktmd2 * phi1max 0 ktmd2];
    Mode = Complex_Eigenvalue_Analysis(M, C, K);
    Damping_ratio2(k2, 1) = min(Mode.("Damping ratio"));
    Frequency2(1:3, k2) = Mode.Frequency;
end

Damping_ratio1_grid = griddata(variables3(:, 1), variables3(:, 2), Damping_ratio1, Xi_syss, Ftmd2s);
Damping_ratio2_grid = griddata(variables3(:, 1), variables3(:, 2), Damping_ratio2, Xi_syss, Ftmd2s);
figure
mesh(Xi_syss, Ftmd2s, Damping_ratio1_grid, 'edgecolor', 'r')
hold on
mesh(Xi_syss, Ftmd2s, Damping_ratio2_grid, 'edgecolor', 'g')

%% 绘制第二个TMD质量比和阻尼比对结构的影响
clc
clear
close all
addpath('..\函数\')
% 结构参数 m c k
m = 1;
k = 27.449810272952913;
omega = sqrt(k / m);
fs = omega / (2 * pi);
% xi_sys = 0.003;
xi_sys=-0.06;
c = 2 * m * omega * xi_sys;
phi1max = 4.765120284718218e-04;

% TMD1参数 mtmd1 ctmd1 ktmd1
mu1 = 0.02;
mtmd1 = mu1 / phi1max^2;
f1 = 1 / (1 + mu1) * fs;
xi1 = sqrt(3 * mu1 / 8 / (1 + mu1));
% f1 = 0.826927;
% xi1 = 0.07;
omega1 = 2 * pi * f1;
ctmd1 = 2 * mtmd1 * omega1 * xi1;
ktmd1 = mtmd1 * omega1^2;

% TMD2参数（变化 质量 和 频率）
mtmd2s = (0.01:0.01:5) * mtmd1;
ftmd2s = (0.7:0.01:1.4) * f1;

% ktmd2s = (0.7:0.01:1.4) * ktmd1;

[Mtmd2s_grid,Ftmd2s_grid]=ndgrid(mtmd2s,ftmd2s);
variables=[Mtmd2s_grid(:) Ftmd2s_grid(:)];

for k1 = 1:size(variables,1)
    xi2=xi1;
    mtmd2=variables(k1,1);
    ftmd2=variables(k1,2);
    omega2=2*pi*ftmd2;
    ktmd2=mtmd2*omega2^2;
    ctmd2=2*mtmd2*omega2*xi2;
    M = [m 0 0; 0 mtmd1 0; 0 0 mtmd2];
    C = [c + ctmd1 * phi1max^2 + ctmd2 * phi1max^2 -ctmd1 * phi1max -ctmd2 * phi1max; -ctmd1 * phi1max ctmd1 0; -ctmd2 * phi1max 0 ctmd2];
    K = [k + ktmd1 * phi1max^2 + ktmd2 * phi1max^2 -ktmd1 * phi1max -ktmd2 * phi1max; -ktmd1 * phi1max ktmd1 0; -ktmd2 * phi1max 0 ktmd2];
    Mode = Complex_Eigenvalue_Analysis(M, C, K);
    Damping_ratio(k1, 1) = min(Mode.("Damping ratio"));
    Frequency(1:3, k1) = Mode.Frequency;
end

Damping_ratio_grid = griddata(variables(:, 1), variables(:, 2), Damping_ratio, Mtmd2s_grid, Ftmd2s_grid);
figure
contour(Mtmd2s_grid, Ftmd2s_grid, Damping_ratio_grid, 'edgecolor', 'r')

testingcases="xi_sys_"+xi_sys;
str="save "+testingcases+".mat "+"xi_sys Mtmd2s_grid Ftmd2s_grid Damping_ratio_grid";
eval(str)



%% 绘制一个TMD质量比和阻尼比对结构的影响
clc
clear
close all
addpath('..\函数\')
% 结构参数 m c k
m = 1;
k = 27.449810272952913;
omega = sqrt(k / m);
fs = omega / (2 * pi);
xi_sys = 0.003;
% xi_sys=-0.005;
c = 2 * m * omega * xi_sys;
phi1max = 4.765120284718218e-04;

% TMD1参数 mtmd1 ctmd1 ktmd1
mu1 = 0.02;
mtmd1 = mu1 / phi1max^2;
f1 = 1 / (1 + mu1) * fs;
xi1 = sqrt(3 * mu1 / 8 / (1 + mu1));
% f1 = 0.826927;
% xi1 = 0.07;
omega1 = 2 * pi * f1;
ctmd1 = 2 * mtmd1 * omega1 * xi1;
ktmd1 = mtmd1 * omega1^2;

% TMD1参数（变化 质量 和 频率）
mtmd1s = (0.01:0.01:5) * mtmd1;
ftmd1s = (0.7:0.01:1.4) * f1;



[Mtmd1s_grid,Ftmd1s_grid]=ndgrid(mtmd1s,ftmd1s);
variables=[Mtmd1s_grid(:) Ftmd1s_grid(:)];

for k1 = 1:size(variables,1)
    mtmd1=variables(k1,1);
    ftmd1=variables(k1,2);
    omega1=2*pi*ftmd1;
    ktmd1=mtmd1*omega1^2;
    ctmd1=2*mtmd1*omega1*xi1;
    M = [m 0; 0 mtmd1;];
    C = [c + ctmd1 * phi1max^2  -ctmd1 * phi1max ; -ctmd1 * phi1max ctmd1];
    K = [k + ktmd1 * phi1max^2  -ktmd1 * phi1max ; -ktmd1 * phi1max ktmd1];
    Mode = Complex_Eigenvalue_Analysis(M, C, K);
    Damping_ratio(k1, 1) = min(Mode.("Damping ratio"));
    Frequency(1:2, k1) = Mode.Frequency;
end

Damping_ratio_grid = griddata(variables(:, 1), variables(:, 2), Damping_ratio, Mtmd1s_grid, Ftmd1s_grid);
figure
contour(Mtmd1s_grid, Ftmd1s_grid, Damping_ratio_grid, 'edgecolor', 'r')

testingcases="single_TMD_xi_sys_"+xi_sys;
str="save "+testingcases+".mat "+"xi_sys Mtmd1s_grid Ftmd1s_grid Damping_ratio_grid";
eval(str)