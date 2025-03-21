%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%Date: 2023-02-03 10:39:23
%LastEditors: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%LastEditTime: 2023-02-06 17:13:36
%FilePath: \NonlinearScanlan\20230124优化问题\粒子群优化\data_analysis.m
%Description: 分析优化过程及结果
%
%Copyright (c) 2023 by ${git_name_email}, All Rights Reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 六个模态六个TMD 最优参数2%质量比
% 分析优化结果
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
mode_re=modeinfo.mode_re;
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
        mode2nodes = mode_re(index(1:2), 1:nModes); %获取两个点坐标的y值
        phi_result = interp1(xResult, mode2nodes, xTMD(t1), 'linear', 'extrap'); %插值以后任意点的振型
        %         disp(phi_result)
        phiTMD_opt(t1, t2) = phi_result(t2);

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
% dis_all_modes=result.dis_all_modes
% save("Dis_with_TMD_opt_177.mat", "dis_all_modes_opt")
% %具体选择哪一行需要修改a_0_main函数中导入的气动力参数 数字表示折减风速（1.77）
% save("Dis_with_TMD_opt_231.mat", "dis_all_modes_opt")
% %具体选择哪一行需要修改a_0_main函数中导入的气动力参数 数字表示折减风速（2.31）
% 生成振型图
figure(2);
plot(nodegap, mode_re(:, 1:6), 'LineWidth', 1.5);
hold on 
phiTMD_plot_opt=[phiTMD_opt(1,1) phiTMD_opt(2,2) phiTMD_opt(3,3) phiTMD_opt(4,4) phiTMD_opt(5,5) phiTMD_opt(6,6)];
scatter(xTMD_opt, phiTMD_plot_opt, 'filled');

%% 六个模态六个TMD 最优参数2%质量比(按照单模态设计)

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

TMDs_location = [276,606,608,394,166,386]; % TMD位置 The location of TMDs
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

mu1=TMDs_mass(1)*phiTMD_original(1)^2;
mu2=TMDs_mass(2)*phiTMD_original(2)^2;
mu3=TMDs_mass(3)*phiTMD_original(3)^2;
mu4=TMDs_mass(4)*phiTMD_original(4)^2;
mu5=TMDs_mass(5)*phiTMD_original(5)^2;
mu6=TMDs_mass(6)*phiTMD_original(6)^2;

f1=1/(1+mu1)*Freq(1);
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
% TMDs_frequency = [0.83, 0.90,1.06,1.28,1.51,1.7]; % TMD频率 The frequency of TMDs
% TMDs_damping_ratio = [0.08, 0.08, 0.08, 0.08, 0.08, 0.08]; % TMD阻尼 The damping of TMDs


% phiTMD 行：TMD的位置 列：TMD位置的模态振型
for t1 = 1:nTMD

    for t2 = 1:nModes
        [~, index] = sort(abs(nodegap - TMDs_location(t1))); %查找与xTMD最接近的点的排序
        xResult = nodegap(index(1:2)); %获取最接近的两个点的x坐标
        mode2nodes = mode_re(index(1:2), 1:nModes); %获取两个点坐标的y值
        phi_result = interp1(xResult, mode2nodes, xTMD(t1), 'linear', 'extrap'); %插值以后任意点的振型
        %         disp(phi_result)
        phiTMD(t1, t2) = phi_result(t2);

    end

end

phiTMD_plot=[phiTMD(1,1) phiTMD(2,2) phiTMD(3,3) phiTMD(4,4) phiTMD(5,5) phiTMD(6,6)];

[result] = a_0_main(number_of_modes_to_control, number_of_modes_to_consider, number_of_tmds, modal_damping_ratios, t_length, TMDs_mass, TMDs_frequency, TMDs_damping_ratio, TMDs_location);
dis_all_modes=result.dis_all_modes
% save("Dis_with_TMD_mode_by_mode_177.mat", "dis_all_modes")
% %具体选择哪一行需要修改a_0_main函数中导入的气动力参数 数字表示折减风速（1.77）
% save("Dis_with_TMD_mode_by_mode_231.mat", "dis_all_modes")
% %具体选择哪一行需要修改a_0_main函数中导入的气动力参数 数字表示折减风速（2.31）

%% 参数对比 Parameters comparison
para_comp=[TMDs_mass' mTMD_opt' TMDs_frequency fTMD_opt' TMDs_damping_ratio' zetaTMD_opt' TMDs_location' xTMD_opt'];
str="质量比为2%TMD，按模态设计以及统一优化的参数绘图信息";
save("6TMD_mu_002_layout.mat","str","para_comp","nodegap","mode_re","phiTMD_plot","phiTMD_plot_opt","TMDs_location","TMDs_damping_ratio","TMDs_frequency","m_oneTMD")
a=1
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


%% 计算振幅折减率（分别为1.77 和 2.31 的折减风速）
clc
clear
close all

Dis_without_TMD_177 = importdata('Dis_without_TMD_177.mat');
Dis_without_TMD_231 = importdata('Dis_without_TMD_231.mat');

Dis_with_TMD_opt_177 = importdata('Dis_with_TMD_opt_177.mat');
Dis_with_TMD_opt_231 = importdata('Dis_with_TMD_opt_231.mat');

Dis_with_TMD_mode_by_mode_177 = importdata('Dis_with_TMD_mode_by_mode_177.mat');
Dis_with_TMD_mode_by_mode_231 = importdata('Dis_with_TMD_mode_by_mode_231.mat');

Dis_without_TMD_177 = Dis_without_TMD_177(1:6);
Dis_without_TMD_231 = Dis_without_TMD_231(1:6);
Dis_with_TMD_opt_177 = Dis_with_TMD_opt_177(1:6);
Dis_with_TMD_opt_231 = Dis_with_TMD_opt_231(1:6);
Dis_with_TMD_mode_by_mode_177 = Dis_with_TMD_mode_by_mode_177(1:6);
Dis_with_TMD_mode_by_mode_231 = Dis_with_TMD_mode_by_mode_231(1:6);

sum_Dis_without_TMD_177 = sum(Dis_without_TMD_177);
sum_Dis_without_TMD_231 = sum(Dis_without_TMD_231);
sum_Dis_with_TMD_opt_177 = sum(Dis_with_TMD_opt_177);
sum_Dis_with_TMD_opt_231 = sum(Dis_with_TMD_opt_231);
sum_Dis_with_TMD_mode_by_mode_177 = sum(Dis_with_TMD_mode_by_mode_177);
sum_Dis_with_TMD_mode_by_mode_231 = sum(Dis_with_TMD_mode_by_mode_231);


modes=[1,2,3,4,5,6];
X = categorical({'Mode 1','Mode 2','Mode 3','Mode 4','Mode 5','Mode 6','Sum'});
modes = reordercats(X,{'Mode 1','Mode 2','Mode 3','Mode 4','Mode 5','Mode 6','Sum'});
x= modes;

y177_mode=1-Dis_with_TMD_mode_by_mode_177./Dis_without_TMD_177;
y177_opt=1-Dis_with_TMD_opt_177./Dis_without_TMD_177;
y231_mode=1-Dis_with_TMD_mode_by_mode_231./Dis_without_TMD_231;
y231_opt=1-Dis_with_TMD_opt_231./Dis_without_TMD_231;

y177_mode=[y177_mode;1-sum_Dis_with_TMD_mode_by_mode_177/sum_Dis_without_TMD_177];
y177_opt=[y177_opt;1-sum_Dis_with_TMD_opt_177/sum_Dis_without_TMD_177];
y231_mode=[y231_mode;1-sum_Dis_with_TMD_mode_by_mode_231/sum_Dis_without_TMD_231];
y231_opt=[y231_opt;1-sum_Dis_with_TMD_opt_231/sum_Dis_without_TMD_231];


figure(1)
bar(x,[y177_mode y177_opt],'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0.5 0.5 0.5])

figure(2)
bar(x,[y231_mode y231_opt],'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0.5 0.5 0.5])
save("Dis_compare_6tmd.mat","x","y177_mode","y177_opt","y231_mode","y231_opt")

%% 六个模态5个TMD 最优参数2%质量比

% 分析优化结果
clc; clear; close all;
addpath('..\')
addpath('..\..\函数\')

% 读取数据
traverse_data = importdata('5TMD粒子群优化\opt_5tmd_1.txt');

% 生成数据表
data_table = array2table(traverse_data, 'VariableNames', {'fTMD1', 'fTMD2', 'fTMD3', 'fTMD4', 'fTMD5', 'mTMD1', 'mTMD2', 'mTMD3', 'mTMD4','xTMD1', 'xTMD2', 'xTMD3', 'xTMD4', 'xTMD5','zetaTMD1', 'zetaTMD2', 'zetaTMD3', 'zetaTMD4', 'zetaTMD5','bestfval','iteration'});

% 生成散点图
figure(1);
scatter(data_table.iteration, data_table.bestfval, 'filled');


Best_Parameters = data_table(data_table.bestfval == min(data_table.bestfval),:);

% 导入桥梁模态信息
modeinfo = load('modeinfo_all.mat'); % 读取未安装TMD时的模态信息. Read the mode information. 
nodegap=modeinfo.nodegap;
mode_re=modeinfo.mode_re;
Freq=modeinfo.Freq;

nTMD = 5;
nModes = 10;
xTMD = Best_Parameters{1, 10:14};
zetaTMD = Best_Parameters{1, 15:19};
fTMD = Best_Parameters{1, 1:5};
mTMD = Best_Parameters{1, 6:9};
total_tmd_mass_ratio = 0.02; % 总质量比 The total mass ratio
mass_six_span = 10007779.7; % 深中通道非通航桥六跨连续梁质量 The mass of 6-span continuous beam of the non-navigational bridge of the Zhenzhong-Link
total_tmd_mass = total_tmd_mass_ratio * mass_six_span; % 总质量 The total mass
mTMD = [mTMD total_tmd_mass-sum(mTMD)]; % 计算最后一个TMD的质量 Calculate the mass of the last TMD

[~, index] = sort(fTMD); % 按照频率从小到大排序 Sort by frequency from small to large
fTMD = fTMD(index);
mTMD = mTMD(index);
xTMD = xTMD(index);
zetaTMD = zetaTMD(index);


para_collect=[mTMD' xTMD' fTMD' zetaTMD'];
save("5TMD_mu_002_layout.mat","para_collect")
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




number_of_modes_to_control=[1 2 3 4 5 6];
number_of_modes_to_consider=10;
number_of_tmds=5;
modal_damping_ratios=ones(1,10)*0.003;
t_length=150;
TMDs_mass=mTMD;
TMDs_frequency=fTMD;
TMDs_damping_ratio=zetaTMD;
TMDs_location=xTMD;


[result] = a_0_main(number_of_modes_to_control,number_of_modes_to_consider,number_of_tmds,modal_damping_ratios,t_length,TMDs_mass,TMDs_frequency,TMDs_damping_ratio,TMDs_location);
dis_all_modes=result.dis_all_modes;

% save("Dis_with_5TMD_opt_177.mat", "dis_all_modes")

% save("Dis_with_5TMD_opt_231.mat", "dis_all_modes")

%% 六个模态3个TMD 最优参数2%质量比

% 分析优化结果
clc; clear; close all;
addpath('..\')
addpath('..\..\函数\')

% 读取数据
traverse_data = importdata('3TMD粒子群优化\opt_3tmd_1.txt');

% 生成数据表
data_table = array2table(traverse_data, 'VariableNames', {'fTMD1', 'fTMD2', 'fTMD3', 'mTMD1', 'mTMD2', 'xTMD1', 'xTMD2', 'xTMD3','zetaTMD1', 'zetaTMD2', 'zetaTMD3','bestfval','iteration'});

% % 生成散点图
% figure(1);
% scatter(data_table.iteration, data_table.bestfval, 'filled');


Best_Parameters = data_table(data_table.bestfval == min(data_table.bestfval),:);

% 导入桥梁模态信息
modeinfo = load('modeinfo_all.mat'); % 读取未安装TMD时的模态信息. Read the mode information. 
nodegap=modeinfo.nodegap;
mode_re=modeinfo.mode_re;
Freq=modeinfo.Freq;

nTMD = 3;
nModes = 10;
xTMD = Best_Parameters{1, 6:8};
zetaTMD = Best_Parameters{1, 9:11};
fTMD = Best_Parameters{1, 1:3};
mTMD = Best_Parameters{1, 4:5};
total_tmd_mass_ratio = 0.02; % 总质量比 The total mass ratio
mass_six_span = 10007779.7; % 深中通道非通航桥六跨连续梁质量 The mass of 6-span continuous beam of the non-navigational bridge of the Zhenzhong-Link
total_tmd_mass = total_tmd_mass_ratio * mass_six_span; % 总质量 The total mass
mTMD = [mTMD total_tmd_mass-sum(mTMD)]; % 计算最后一个TMD的质量 Calculate the mass of the last TMD

[~, index] = sort(fTMD); % 按照频率从小到大排序 Sort by frequency from small to large
fTMD = fTMD(index);
mTMD = mTMD(index);
xTMD = xTMD(index);
zetaTMD = zetaTMD(index);


para_collect=[mTMD' xTMD' fTMD' zetaTMD'];
save("3TMD_mu_002_layout.mat","para_collect")


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




number_of_modes_to_control=[1 2 3 4 5 6];
number_of_modes_to_consider=10;
number_of_tmds=3;
modal_damping_ratios=ones(1,10)*0.003;
t_length=150;
TMDs_mass=mTMD;
TMDs_frequency=fTMD;
TMDs_damping_ratio=zetaTMD;
TMDs_location=xTMD;

[result] = a_0_main(number_of_modes_to_control,number_of_modes_to_consider,number_of_tmds,modal_damping_ratios,t_length,TMDs_mass,TMDs_frequency,TMDs_damping_ratio,TMDs_location);
dis_all_modes=result.dis_all_modes;

% save("Dis_with_3TMD_opt_177.mat", "dis_all_modes")

% save("Dis_with_3TMD_opt_231.mat", "dis_all_modes")


%% 计算振幅折减率（分别为1.77 和 2.31 的折减风速）
clc
clear
close all

Dis_without_TMD_177 = importdata('Dis_without_TMD_177.mat');
Dis_without_TMD_231 = importdata('Dis_without_TMD_231.mat');

Dis_with_TMD_opt_177 = importdata('Dis_with_TMD_opt_177.mat');
Dis_with_TMD_opt_231 = importdata('Dis_with_TMD_opt_231.mat');
Dis_with_5TMD_opt_177 = importdata('Dis_with_5TMD_opt_177.mat');
Dis_with_5TMD_opt_231 = importdata('Dis_with_5TMD_opt_231.mat');
Dis_with_3TMD_opt_177 = importdata('Dis_with_3TMD_opt_177.mat');
Dis_with_3TMD_opt_231 = importdata('Dis_with_3TMD_opt_231.mat');


Dis_with_TMD_mode_by_mode_177 = importdata('Dis_with_TMD_mode_by_mode_177.mat');
Dis_with_TMD_mode_by_mode_231 = importdata('Dis_with_TMD_mode_by_mode_231.mat');

Dis_without_TMD_177 = Dis_without_TMD_177(1:6);
Dis_without_TMD_231 = Dis_without_TMD_231(1:6);
Dis_with_TMD_opt_177 = Dis_with_TMD_opt_177(1:6);
Dis_with_TMD_opt_231 = Dis_with_TMD_opt_231(1:6);
Dis_with_TMD_mode_by_mode_177 = Dis_with_TMD_mode_by_mode_177(1:6);
Dis_with_TMD_mode_by_mode_231 = Dis_with_TMD_mode_by_mode_231(1:6);
Dis_with_5TMD_opt_177=Dis_with_5TMD_opt_177(1:6);
Dis_with_5TMD_opt_231=Dis_with_5TMD_opt_231(1:6);
Dis_with_3TMD_opt_177=Dis_with_3TMD_opt_177(1:6);
Dis_with_3TMD_opt_231=Dis_with_3TMD_opt_231(1:6);

sum_Dis_without_TMD_177 = sum(Dis_without_TMD_177);
sum_Dis_without_TMD_231 = sum(Dis_without_TMD_231);
sum_Dis_with_TMD_opt_177 = sum(Dis_with_TMD_opt_177);
sum_Dis_with_TMD_opt_231 = sum(Dis_with_TMD_opt_231);
sum_Dis_with_TMD_mode_by_mode_177 = sum(Dis_with_TMD_mode_by_mode_177);
sum_Dis_with_TMD_mode_by_mode_231 = sum(Dis_with_TMD_mode_by_mode_231);
sum_Dis_with_5TMD_opt_177=sum(Dis_with_5TMD_opt_177);
sum_Dis_with_5TMD_opt_231=sum(Dis_with_5TMD_opt_231);
sum_Dis_with_3TMD_opt_177=sum(Dis_with_3TMD_opt_177);
sum_Dis_with_3TMD_opt_231=sum(Dis_with_3TMD_opt_231);



modes=[1,2,3,4,5,6];
X = categorical({'Mode 1','Mode 2','Mode 3','Mode 4','Mode 5','Mode 6','Sum'});
modes = reordercats(X,{'Mode 1','Mode 2','Mode 3','Mode 4','Mode 5','Mode 6','Sum'});
x= modes;

y177_mode=1-Dis_with_TMD_mode_by_mode_177./Dis_without_TMD_177;
y177_opt=1-Dis_with_TMD_opt_177./Dis_without_TMD_177;
y231_mode=1-Dis_with_TMD_mode_by_mode_231./Dis_without_TMD_231;
y231_opt=1-Dis_with_TMD_opt_231./Dis_without_TMD_231;

y177_opt_5TMD=1-Dis_with_5TMD_opt_177./Dis_without_TMD_177;
y231_opt_5TMD=1-Dis_with_5TMD_opt_231./Dis_without_TMD_231;
y177_opt_3TMD=1-Dis_with_3TMD_opt_177./Dis_without_TMD_177;
y231_opt_3TMD=1-Dis_with_3TMD_opt_231./Dis_without_TMD_231;

y177_mode=[y177_mode;1-sum_Dis_with_TMD_mode_by_mode_177/sum_Dis_without_TMD_177];
y177_opt=[y177_opt;1-sum_Dis_with_TMD_opt_177/sum_Dis_without_TMD_177];
y231_mode=[y231_mode;1-sum_Dis_with_TMD_mode_by_mode_231/sum_Dis_without_TMD_231];
y231_opt=[y231_opt;1-sum_Dis_with_TMD_opt_231/sum_Dis_without_TMD_231];

y177_opt_5TMD=[y177_opt_5TMD;1-sum_Dis_with_5TMD_opt_177/sum_Dis_without_TMD_177];
y231_opt_5TMD=[y231_opt_5TMD;1-sum_Dis_with_5TMD_opt_231/sum_Dis_without_TMD_231];
y177_opt_3TMD=[y177_opt_3TMD;1-sum_Dis_with_3TMD_opt_177/sum_Dis_without_TMD_177];
y231_opt_3TMD=[y231_opt_3TMD;1-sum_Dis_with_3TMD_opt_231/sum_Dis_without_TMD_231];


figure(1)
bar(x,[y177_mode y177_opt y177_opt_5TMD y177_opt_3TMD],'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0.5 0.5 0.5])

figure(2)
bar(x,[y231_mode y231_opt y231_opt_5TMD y231_opt_3TMD],'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0.5 0.5 0.5])
save("Dis_compare_tmds.mat","x","y177_mode","y177_opt","y231_mode","y231_opt","y177_opt_5TMD","y177_opt_3TMD","y231_opt_5TMD","y231_opt_3TMD")