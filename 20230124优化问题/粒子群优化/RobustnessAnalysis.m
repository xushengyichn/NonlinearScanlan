clc;clear;close all
addpath('..\')
addpath('..\..\函数\')

%% Scheme 1
six_TMD_mu_002_layout = importdata('6TMD_mu_002_layout.mat');


modeinfo = load('modeinfo_all.mat'); % 读取未安装TMD时的模态信息. Read the mode information. 
nodegap=modeinfo.nodegap;
mode_re=modeinfo.mode_re;

number_of_modes_to_control = [1, 2, 3, 4, 5, 6]; % 需要控制的模态数 The number of modes to be controlled
number_of_modes_to_consider = 10; % 考虑的总模态数 The total number of modes considered
number_of_tmds = 6; % 考虑的总TMD数 The total number of TMDs
nTMD = number_of_tmds;
nModes = number_of_modes_to_consider;
modal_damping_ratios = ones(1, number_of_modes_to_consider) * 0.003; % 模态阻尼 The damping of modes
t_length = 150; % 设定计算时间长度
para_comp = six_TMD_mu_002_layout.para_comp;
TMDs_mass = para_comp(:,1);
TMDs_frequency = para_comp(:,3);
TMDs_damping_ratio = para_comp(:,5);
TMDs_location = para_comp(:,7);
xTMD = TMDs_location;
freqratios = -20:1:20;

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
    S1_data_all(:,k1)=dis_all_modes(1:6,:);
end

figure 
plot(freqratios,S1_data)

%% Scheme 2
six_TMD_mu_002_layout = importdata('6TMD_mu_002_layout.mat');


modeinfo = load('modeinfo_all.mat'); % 读取未安装TMD时的模态信息. Read the mode information. 
nodegap=modeinfo.nodegap;
mode_re=modeinfo.mode_re;

number_of_modes_to_control = [1, 2, 3, 4, 5, 6]; % 需要控制的模态数 The number of modes to be controlled
number_of_modes_to_consider = 10; % 考虑的总模态数 The total number of modes considered
number_of_tmds = 6; % 考虑的总TMD数 The total number of TMDs
nTMD = number_of_tmds;
nModes = number_of_modes_to_consider;
modal_damping_ratios = ones(1, number_of_modes_to_consider) * 0.003; % 模态阻尼 The damping of modes
t_length = 150; % 设定计算时间长度
para_comp = six_TMD_mu_002_layout.para_comp;
TMDs_mass = para_comp(:,2);
TMDs_frequency = para_comp(:,4);
TMDs_damping_ratio = para_comp(:,6);
TMDs_location = para_comp(:,8);
xTMD = TMDs_location;
freqratios = -20:1:20;

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

for k1 = 1:length(freqratios)
    freqratio=freqratios(k1);
    % freqratio=0
    [result] = a_0_robust(number_of_modes_to_control, number_of_modes_to_consider, number_of_tmds, modal_damping_ratios, t_length, TMDs_mass, TMDs_frequency, TMDs_damping_ratio, TMDs_location,freqratio);
    S2_data(k1)=result.dis_all_modes_sum;
    dis_all_modes=result.dis_all_modes;
    S2_data_all(:,k1)=dis_all_modes(1:6,:);
end

figure 
plot(freqratios,S2_data)


%% Scheme 3
five_TMD_mu_002_layout = importdata('5TMD_mu_002_layout.mat');

modeinfo = load('modeinfo_all.mat'); % 读取未安装TMD时的模态信息. Read the mode information. 
nodegap=modeinfo.nodegap;
mode_re=modeinfo.mode_re;

number_of_modes_to_control = [1, 2, 3, 4, 5, 6]; % 需要控制的模态数 The number of modes to be controlled
number_of_modes_to_consider = 10; % 考虑的总模态数 The total number of modes considered
number_of_tmds = 5; % 考虑的总TMD数 The total number of TMDs
nTMD = number_of_tmds;
nModes = number_of_modes_to_consider;
modal_damping_ratios = ones(1, number_of_modes_to_consider) * 0.003; % 模态阻尼 The damping of modes
t_length = 150; % 设定计算时间长度

TMDs_mass = five_TMD_mu_002_layout(:,1);
TMDs_frequency = five_TMD_mu_002_layout(:,3);
TMDs_damping_ratio = five_TMD_mu_002_layout(:,4);
TMDs_location = five_TMD_mu_002_layout(:,2);
xTMD = TMDs_location;
freqratios = -20:1:20;

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

for k1 = 1:length(freqratios)
    freqratio=freqratios(k1);
    % freqratio=0
    [result] = a_0_robust(number_of_modes_to_control, number_of_modes_to_consider, number_of_tmds, modal_damping_ratios, t_length, TMDs_mass, TMDs_frequency, TMDs_damping_ratio, TMDs_location,freqratio);
    S3_data(k1)=result.dis_all_modes_sum;
    dis_all_modes=result.dis_all_modes;
    S3_data_all(:,k1)=dis_all_modes(1:6,:);
end

figure 
plot(freqratios,S3_data)

%% Scheme 4
three_TMD_mu_002_layout = importdata('3TMD_mu_002_layout.mat');

modeinfo = load('modeinfo_all.mat'); % 读取未安装TMD时的模态信息. Read the mode information. 
nodegap=modeinfo.nodegap;
mode_re=modeinfo.mode_re;

number_of_modes_to_control = [1, 2, 3, 4, 5, 6]; % 需要控制的模态数 The number of modes to be controlled
number_of_modes_to_consider = 10; % 考虑的总模态数 The total number of modes considered
number_of_tmds = 3; % 考虑的总TMD数 The total number of TMDs
nTMD = number_of_tmds;
nModes = number_of_modes_to_consider;
modal_damping_ratios = ones(1, number_of_modes_to_consider) * 0.003; % 模态阻尼 The damping of modes
t_length = 150; % 设定计算时间长度

TMDs_mass = three_TMD_mu_002_layout(:,1);
TMDs_frequency = three_TMD_mu_002_layout(:,3);
TMDs_damping_ratio = three_TMD_mu_002_layout(:,4);
TMDs_location = three_TMD_mu_002_layout(:,2);
xTMD = TMDs_location;
freqratios = -20:1:20;

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

for k1 = 1:length(freqratios)
    freqratio=freqratios(k1);
    % freqratio=0
    [result] = a_0_robust(number_of_modes_to_control, number_of_modes_to_consider, number_of_tmds, modal_damping_ratios, t_length, TMDs_mass, TMDs_frequency, TMDs_damping_ratio, TMDs_location,freqratio);
    S4_data(k1)=result.dis_all_modes_sum;
    dis_all_modes=result.dis_all_modes;
    S4_data_all(:,k1)=dis_all_modes(1:6,:);
end

figure 
plot(freqratios,S4_data)

figure
plot(freqratios,S1_data)
hold on
plot(freqratios,S2_data)
plot(freqratios,S3_data)
plot(freqratios,S4_data)

legend('s1','s2','s3','s4')

save("RobustnessAnalysis_freq.mat","freqratios","S1_data","S1_data_all","S2_data","S2_data_all","S3_data","S3_data_all","S4_data","S4_data_all")