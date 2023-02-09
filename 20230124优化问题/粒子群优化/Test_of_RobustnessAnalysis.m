%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%Date: 2023-02-07 13:16:41
%LastEditors: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%LastEditTime: 2023-02-07 13:18:59
%FilePath: \NonlinearScanlan\20230124优化问题\粒子群优化\Test_of_RobustnessAnalysis.m
%Description: 
%
%Copyright (c) 2023 by ${git_name_email}, All Rights Reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear;close all
addpath('..\')
addpath('..\..\函数\')

%% Scheme 1
six_TMD_mu_002_layout = importdata('6TMD_mu_002_layout.mat');


modeinfo = load('modeinfo_all.mat'); % 读取未安装TMD时的模态信息. Read the mode information. 
nodegap=modeinfo.nodegap;
mode_re=modeinfo.mode_re;
mode = modeinfo.mode;
Freq= modeinfo.Freq;

number_of_modes_to_control = [1 2 3 4 5 6]; % 需要控制的模态数 The number of modes to be controlled
number_of_modes_to_consider = 10; % 考虑的总模态数 The total number of modes considered
number_of_tmds = 6; % 考虑的总TMD数 The total number of TMDs
nTMD = number_of_tmds;
nModes = number_of_modes_to_consider;
modal_damping_ratios = ones(1, number_of_modes_to_consider) * 0.03; % 模态阻尼 The damping of modes
t_length = 150; % 设定计算时间长度
para_comp = six_TMD_mu_002_layout.para_comp;
TMDs_mass = para_comp(:,1);

TMDs_location = para_comp(:,7);
xTMD = TMDs_location;

for t1 = 1:nTMD

    for t2 = 1:nModes
        [~, index] = sort(abs(nodegap - xTMD(t1))); %查找与xTMD最接近的点的排序
        xResult = nodegap(index(1:2)); %获取最接近的两个点的x坐标
        mode2nodes = mode(index(1:2), 1:nModes); %获取两个点坐标的y值
        phi_result = interp1(xResult, mode2nodes, xTMD(t1), 'linear', 'extrap'); %插值以后任意点的振型
        %         disp(phi_result)
        phiTMD_original(t1, t2) = phi_result(t2);

    end

end

maxmode= max(mode)
if nTMD >0
mu1=TMDs_mass(1)*phiTMD_original(1,1)^2;
else
    mu1=0
end
mu2=TMDs_mass(2)*phiTMD_original(2,2)^2;
mu3=TMDs_mass(3)*phiTMD_original(3,3)^2;
mu4=TMDs_mass(4)*phiTMD_original(4,4)^2;
mu5=TMDs_mass(5)*phiTMD_original(5,5)^2;
mu6=TMDs_mass(6)*phiTMD_original(6,6)^2;

f1=1/(1+mu1)*Freq(1);
% f1=Freq(1)*1.01;
f2=1/(1+mu2)*Freq(2);
f3=1/(1+mu3)*Freq(3);
f4=1/(1+mu4)*Freq(4);
f5=1/(1+mu5)*Freq(5);
f6=1/(1+mu6)*Freq(6);

zeta1 = sqrt(3*mu1/8/(1+mu1));
zeta2 = sqrt(3*mu2/8/(1+mu2));
zeta3 = sqrt(3*mu3/8/(1+mu3));
zeta4 = sqrt(3*mu4/8/(1+mu4));
zeta5 = sqrt(3*mu5/8/(1+mu5));
zeta6 = sqrt(3*mu6/8/(1+mu6));

TMDs_frequency = [f1 f2 f3 f4 f5 f6]';
TMDs_damping_ratio = [zeta1 zeta2 zeta3 zeta4 zeta5 zeta6];

% TMDs_frequency = [f1]';
% TMDs_damping_ratio = [zeta1];
% TMDs_frequency = para_comp(:,3);
% TMDs_damping_ratio = para_comp(:,5);

freqratios = -5:1:5;

% phiTMD 行：TMD的位置 列：TMD位置的模态振型
for t1 = 1:nTMD

    for t2 = 1:nModes
        [~, index] = sort(abs(nodegap - xTMD(t1))); %查找与xTMD最接近的点的排序
        xResult = nodegap(index(1:2)); %获取最接近的两个点的x坐标
        mode2nodes = mode_re(index(1:2), 1:nModes); %获取两个点坐标的y值
        phi_result = interp1(xResult, mode2nodes, xTMD(t1), 'linear', 'extrap'); %插值以后任意点的振型
        %         disp(phi_result)
        phiTMD_opt(t1, t2) = phi_result(t2);

    end

end

parfor k1 = 1:length(freqratios)
    freqratio=freqratios(k1);
    % freqratio=0
    [result] = a_0_robust(number_of_modes_to_control, number_of_modes_to_consider, number_of_tmds, modal_damping_ratios, t_length, TMDs_mass, TMDs_frequency, TMDs_damping_ratio, TMDs_location,freqratio);
    S1_data(k1)=result.dis_all_modes_sum;
    dis_all_modes=result.dis_all_modes;
    S1_data_all(:,k1)=dis_all_modes(1:1,:);
end

figure 
plot(freqratios,S1_data)

