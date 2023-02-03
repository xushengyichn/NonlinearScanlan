%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%Date: 2023-02-03 10:39:23
%LastEditors: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%LastEditTime: 2023-02-03 13:55:07
%FilePath: \粒子群优化\data_analysis.m
%Description: 分析优化过程及结果
%
%Copyright (c) 2023 by ${git_name_email}, All Rights Reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 分析穷举结果
clc; clear; close all;
addpath('..\')
addpath('..\..\函数\')

% 读取数据
traverse_data = importdata('trial_points_worker_1.txt');

% 生成数据表
data_table = array2table(traverse_data, 'VariableNames', {'fTMD1', 'fTMD2', 'fTMD3', 'fTMD4', 'fTMD5', 'fTMD6', 'mTMD1', 'mTMD2', 'mTMD3', 'mTMD4', 'mTMD5','xTMD1', 'xTMD2', 'xTMD3', 'xTMD4', 'xTMD5', 'xTMD6','zetaTMD1', 'zetaTMD2', 'zetaTMD3', 'zetaTMD4', 'zetaTMD5', 'zetaTMD6','bestfval','iteration'});

% 生成散点图
figure(1);
scatter(data_table.iteration, data_table.bestfval, 'filled');


Best_Parameters = data_table(data_table.bestfval == min(data_table.bestfval),:);

% 导入桥梁模态信息
modeinfo = load('modeinfo_all.mat'); % 读取未安装TMD时的模态信息. Read the mode information. 
nodegap=modeinfo.nodegap;
mode=modeinfo.mode_re;
Freq=modeinfo.Freq;

nTMD = 6;
nModes = 10;
xTMD = Best_Parameters{1, 12:17};
zetaTMD = Best_Parameters{1, 18:23};
fTMD = Best_Parameters{1, 1:6};
mTMD = Best_Parameters{1, 7:11};
total_tmd_mass_ratio = 0.02; % 总质量比 The total mass ratio
mass_six_span = 10007779.7; % 深中通道非通航桥六跨连续梁质量 The mass of 6-span continuous beam of the non-navigational bridge of the Zhenzhong-Link
total_tmd_mass = total_tmd_mass_ratio * mass_six_span; % 总质量 The total mass
mTMD = [mTMD total_tmd_mass-sum(mTMD)]; % 计算最后一个TMD的质量 Calculate the mass of the last TMD

[~, index] = sort(fTMD); % 按照频率从小到大排序 Sort by frequency from small to large
fTMD = fTMD(index);
mTMD = mTMD(index);
xTMD = xTMD(index);
zetaTMD = zetaTMD(index);

% phiTMD 行：TMD的位置 列：TMD位置的模态振型
for t1 = 1:nTMD

    for t2 = 1:nModes
        [~, index] = sort(abs(nodegap - xTMD(t1))); %查找与xTMD最接近的点的排序
        xResult = nodegap(index(1:2)); %获取最接近的两个点的x坐标
        mode2nodes = mode(index(1:2), 1:nModes); %获取两个点坐标的y值
        phi_result = interp1(xResult, mode2nodes, xTMD(t1), 'linear', 'extrap'); %插值以后任意点的振型
        %         disp(phi_result)
        phiTMD(t1, t2) = phi_result(t2);

    end

end




number_of_modes_to_control=[1 2 3 4 5 6];
number_of_modes_to_consider=10;
number_of_tmds=6;
modal_damping_ratios=ones(1,10)*0.003;
t_length=150;
TMDs_mass=mTMD;
TMDs_frequency=fTMD;
TMDs_damping_ratio=zetaTMD;
TMDs_location=xTMD;
[result] = a_0_main(number_of_modes_to_control,number_of_modes_to_consider,number_of_tmds,modal_damping_ratios,t_length,TMDs_mass,TMDs_frequency,TMDs_damping_ratio,TMDs_location);
dis_all_modes_opt=result.dis_all_modes;
xTMD_opt =xTMD;
zetaTMD_opt =zetaTMD;
fTMD_opt =fTMD;
mTMD_opt =mTMD;

% 生成振型图
figure(2);
plot(nodegap, mode(:, 1:6), 'LineWidth', 1.5);
hold on 
phiTMD_plot_opt=[phiTMD(1,1) phiTMD(2,2) phiTMD(3,3) phiTMD(4,4) phiTMD(5,5) phiTMD(6,6)];
scatter(xTMD_opt, phiTMD_plot_opt, 'filled');

%% 六个模态六个TMD 最优参数2%质量比

% Preset parameters

number_of_modes_to_control = [1, 2, 3, 4, 5, 6]; % 需要控制的模态数 The number of modes to be controlled
number_of_modes_to_consider = 10; % 考虑的总模态数 The total number of modes considered
number_of_tmds = 6; % 考虑的总TMD数 The total number of TMDs
total_tmd_mass_ratio = 0.02; % 总质量比 The total mass ratio
mass_six_span = 10007779.7; % 深中通道非通航桥六跨连续梁质量 The mass of 6-span continuous beam of the non-navigational bridge of the Zhenzhong-Link
total_tmd_mass = total_tmd_mass_ratio * mass_six_span; % 总质量 The total mass
modal_damping_ratios = ones(1, number_of_modes_to_consider) * 0.003; % 模态阻尼 The damping of modes
t_length = 100; % 设定计算时间长度

% 输入初始参数
% Input initial parameters
m_oneTMD= total_tmd_mass / number_of_tmds;
TMDs_mass = [m_oneTMD,m_oneTMD,m_oneTMD,m_oneTMD,m_oneTMD]; % TMD质量 The mass of TMDs
TMDs_mass = [TMDs_mass, total_tmd_mass - sum(TMDs_mass)]; % 最后一个TMD的质量 The mass of the last TMD

% 计算模态质量
modeinfo = load('modeinfo_all.mat'); % 读取未安装TMD时的模态信息. Read the mode information. 
nodegap=modeinfo.nodegap;
mode=modeinfo.mode;
[modemax, modemaxloc]=max(abs(mode(:,1:6)));
mu=m_oneTMD.*modemax.^2;
locations=nodegap(modemaxloc);


TMDs_frequency = [0.83, 0.90,1.06,1.28,1.51,1.7]; % TMD频率 The frequency of TMDs
TMDs_damping_ratio = [0.08, 0.08, 0.08, 0.08, 0.08, 0.08]; % TMD阻尼 The damping of TMDs
TMDs_location = [276,606,608,394,166,386]; % TMD位置 The location of TMDs

[result] = a_0_main(number_of_modes_to_control, number_of_modes_to_consider, number_of_tmds, modal_damping_ratios, t_length, TMDs_mass, TMDs_frequency, TMDs_damping_ratio, TMDs_location);
dis_all_modes=result.dis_all_modes


%% 参数对比 Parameters comparison
para_comp=[TMDs_mass' mTMD_opt' TMDs_frequency' fTMD_opt' TMDs_damping_ratio' zetaTMD_opt' TMDs_location' xTMD_opt'];


%% 不安装TMD的原始振幅

% Preset parameters

number_of_modes_to_control = [1, 2, 3, 4, 5, 6]; % 需要控制的模态数 The number of modes to be controlled
number_of_modes_to_consider = 10; % 考虑的总模态数 The total number of modes considered
number_of_tmds = 0; % 考虑的总TMD数 The total number of TMDs
total_tmd_mass_ratio = 0.02; % 总质量比 The total mass ratio
mass_six_span = 10007779.7; % 深中通道非通航桥六跨连续梁质量 The mass of 6-span continuous beam of the non-navigational bridge of the Zhenzhong-Link
total_tmd_mass = total_tmd_mass_ratio * mass_six_span; % 总质量 The total mass
modal_damping_ratios = ones(1, number_of_modes_to_consider) * 0.003; % 模态阻尼 The damping of modes
t_length = 100; % 设定计算时间长度

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

[result] = a_0_main(number_of_modes_to_control, number_of_modes_to_consider, number_of_tmds, modal_damping_ratios, t_length, TMDs_mass, TMDs_frequency, TMDs_damping_ratio, TMDs_location);
dis_all_modes=result.dis_all_modes

% save("Dis_without_TMD_177.m","dis_all_modes")
% %具体选择哪一行需要修改a_0_main函数中导入的气动力参数 数字表示折减风速（1.77）
% save("Dis_without_TMD_231.m","dis_all_modes")
% %具体选择哪一行需要修改a_0_main函数中导入的气动力参数 数字表示折减风速（2.31）