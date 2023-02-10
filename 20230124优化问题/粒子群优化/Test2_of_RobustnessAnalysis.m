clc;clear;close all
addpath('..\')
addpath('..\..\函数\')


damping_ratio_ratios = -5;

%% Scheme 4






three_TMD_mu_002_layout = importdata('3TMD_mu_002_layout.mat');

modeinfo = load('modeinfo_all.mat'); % 读取未安装TMD时的模态信息. Read the mode information. 
nodegap=modeinfo.nodegap;
mode_re=modeinfo.mode_re;
mode = modeinfo.mode;

[phideckmax,location] = max(mode); %计算桥面每个模态的最大振型值

eig_val = modeinfo.eig_val;
eig_vec = modeinfo.eig_vec;
[omeg, w_order] = sort(sqrt(diag(eig_val)));
mode_vec = eig_vec(:, w_order);
Freq = omeg / (2 * pi);


for k1 = 1:length(nodegap) - 1
    nodedis(k1, 1) = nodegap(k1 + 1) - nodegap(k1);
end

% 模态振型比节点间距多一个点，所以将两节点的模态取平均值计算振型积分
for k1 = 1:length(nodegap) - 1
    modecal(k1, :) = (mode(k1 + 1, :) + mode(k1, :)) / 2;
end




number_of_modes_to_control = [1]; % 需要控制的模态数 The number of modes to be controlled
number_of_modes_to_consider = 10; % 考虑的总模态数 The total number of modes considered
number_of_tmds = 3; % 考虑的总TMD数 The total number of TMDs
nTMD = number_of_tmds;
nModes = number_of_modes_to_consider;
modal_damping_ratios = ones(1, number_of_modes_to_consider) * 0.003; % 模态阻尼 The damping of modes
t_length = 1000; % 设定计算时间长度




for k1 = 1:nModes
    integral_1(k1) = sum(abs(modecal(:, k1)) .* nodedis);
    integral_2(k1) = sum(modecal(:, k1).^2 .* nodedis);
    integral_3(k1) = sum(modecal(:, k1).^2 .* abs(modecal(:, k1)) .* nodedis);
    integral_4(k1) = sum(modecal(:, k1).^4 .* nodedis);
    integral_5(k1) = sum(modecal(:, k1).^2 .* abs(modecal(:, k1)).^3 .* nodedis);
    integral_6(k1) = sum(modecal(:, k1).^6 .* nodedis);
end

mode_number = number_of_modes_to_control;
phi = mode(:, mode_number);
mode_integral_1 = integral_1(mode_number);
mode_integral_2 = integral_2(mode_number);
mode_integral_3 = integral_3(mode_number);
mode_integral_4 = integral_4(mode_number);
mode_integral_5 = integral_5(mode_number);
mode_integral_6 = integral_6(mode_number);


TMDs_mass = three_TMD_mu_002_layout(:,1);
TMDs_frequency = three_TMD_mu_002_layout(:,3);
TMDs_damping_ratio = three_TMD_mu_002_layout(:,4);
TMDs_location = three_TMD_mu_002_layout(:,2);
xTMD = TMDs_location;

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

for k1 = 1:length(damping_ratio_ratios)
    damping_ratio_ratios_temp=damping_ratio_ratios(k1);
    modal_damping_ratios=ones(1, number_of_modes_to_consider) * 0.003*(1+damping_ratio_ratios_temp/100)
    freqratio=0
    [result] = a_0_robust(number_of_modes_to_control, number_of_modes_to_consider, number_of_tmds, modal_damping_ratios, t_length, TMDs_mass, TMDs_frequency, TMDs_damping_ratio, TMDs_location,freqratio);
    S4_data(k1)=result.dis_all_modes_sum;
    dis_all_modes=result.dis_all_modes;
    S4_data_all(:,k1)=dis_all_modes(1:6,:);
end
dis_TMD=result.dis_TMD;

t=result.t;
seqs_allmodes_max_point=result.seqs_allmodes_max_point;

subplot(2,2,1)
plot(t,seqs_allmodes_max_point)

subplot(2,2,2)
plot(t,dis_TMD(:,1))

subplot(2,2,3)
plot(t,dis_TMD(:,2))

subplot(2,2,4)
plot(t,dis_TMD(:,3))

h=t(2)-t(1);
fs=1/h;
time_history_cut=seqs_allmodes_max_point(round(end / 8 * 6, 0):end);

[psd_avg, f, psd_plot] = fft_transfer(fs,time_history_cut');
figure
plot(f,psd_plot)

MM = result.MM;
CC = result.CC;
KK = result.KK;

CC_aero=CC;
rho=result.rho;
U=result.U;
a1=result.a1;
D=result.D;
CC_aero(mode_number,mode_number)=CC_aero(mode_number,mode_number)-rho*U*D*a1*mode_integral_2;

Mode = Complex_Eigenvalue_Analysis(MM,CC_aero,KK);

% 
% figure 
% plot(damping_ratio_ratios,S4_data)

% figure
% plot(damping_ratio_ratios,S1_data)
% hold on
% plot(damping_ratio_ratios,S2_data)
% plot(damping_ratio_ratios,S3_data)
% plot(damping_ratio_ratios,S4_data)

legend('s1','s2','s3','s4')

% save("RobustnessAnalysis_damp.mat","damping_ratio_ratios","S1_data","S1_data_all","S2_data","S2_data_all","S3_data","S3_data_all","S4_data","S4_data_all")