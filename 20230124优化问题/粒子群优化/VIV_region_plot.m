%% 初始化
clc; clear; close all;
addpath('..\')
addpath('..\..\函数\')



%% 不安装TMD的原始振幅

% Preset parameters

number_of_modes_to_control = [1, 2, 3, 4, 5, 6]; % 需要控制的模态数 The number of modes to be controlled
number_of_modes_to_consider = 6; % 考虑的总模态数 The total number of modes considered
number_of_tmds = 0; % 考虑的总TMD数 The total number of TMDs
total_tmd_mass_ratio = 0.02; % 总质量比 The total mass ratio
mass_six_span = 10007779.7; % 深中通道非通航桥六跨连续梁质量 The mass of 6-span continuous beam of the non-navigational bridge of the Zhenzhong-Link
total_tmd_mass = total_tmd_mass_ratio * mass_six_span; % 总质量 The total mass
modal_damping_ratios = ones(1, number_of_modes_to_consider) * 0.003; % 模态阻尼 The damping of modes
t_length = 200; % 设定计算时间长度

% 输入初始参数
% Input initial parameters
m_oneTMD= total_tmd_mass / number_of_tmds;
TMDs_mass = [m_oneTMD,m_oneTMD,m_oneTMD,m_oneTMD,m_oneTMD]*0; % TMD质量 The mass of TMDs
TMDs_mass = [TMDs_mass, total_tmd_mass - sum(TMDs_mass)]*0; % 最后一个TMD的质量 The mass of the last TMD

% 计算模态质量
modeinfo = load('modeinfo_all.mat'); % 读取未安装TMD时的模态信息. Read the mode information. 
nodegap=modeinfo.nodegap;
mode=modeinfo.mode;
[modemax, modemaxloc]=max(abs(mode(:,1:6)));
mu=m_oneTMD.*modemax.^2;
locations=nodegap(modemaxloc);


TMDs_frequency = [0.83, 0.90,1.06,1.28,1.51,1.7]*0; % TMD频率 The frequency of TMDs
TMDs_damping_ratio = [0.08, 0.08, 0.08, 0.08, 0.08, 0.08]*0; % TMD阻尼 The damping of TMDs
TMDs_location = [276,606,608,394,166,386]*0; % TMD位置 The location of TMDs

VIV_windspeeds = [1,2,3,5,6,7,8,10];

for k1 = 1:length(number_of_modes_to_control)
    for k2 =1:length(VIV_windspeeds)
        VIV_windspeed=VIV_windspeeds(k2);
        number_of_modes_to_control_temp=number_of_modes_to_control(k1);
        [result] = d_0_TMD_without_TMD(number_of_modes_to_control_temp, number_of_modes_to_consider, number_of_tmds, modal_damping_ratios, t_length, TMDs_mass, TMDs_frequency, TMDs_damping_ratio, TMDs_location,VIV_windspeed);
        dis_all_modes=result.dis_all_modes;
        dis_all_modes_sum(k1,k2) = result.dis_all_modes_sum;
        U(k1,k2)=result.U;
        a=1;
    end
end

figure
hold on
for k1 = 1:6
    plot(U(k1,:),dis_all_modes_sum(k1,:))
end

save("VIV_region_plot.mat","U","dis_all_modes_sum")