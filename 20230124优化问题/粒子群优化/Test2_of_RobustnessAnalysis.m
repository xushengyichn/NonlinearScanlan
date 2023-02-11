clc;clear;close all
addpath('..\')
addpath('..\..\函数\')


damping_ratio_ratios = 32;

%% Scheme 2
six_TMD_mu_002_layout = importdata('6TMD_mu_002_layout.mat');


modeinfo = load('modeinfo_all.mat'); % 读取未安装TMD时的模态信息. Read the mode information. 
nodegap=modeinfo.nodegap;
mode_re=modeinfo.mode_re;

% number_of_modes_to_control = [1, 2, 3, 4, 5, 6]; % 需要控制的模态数 The number of modes to be controlled
number_of_modes_to_control = [1]; % 需要控制的模态数 The number of modes to be controlled
number_of_modes_to_consider = 14; % 考虑的总模态数 The total number of modes considered
number_of_tmds = 6; % 考虑的总TMD数 The total number of TMDs
nTMD = number_of_tmds;
nModes = number_of_modes_to_consider;
modal_damping_ratios = ones(1, number_of_modes_to_consider) * 0.003; % 模态阻尼 The damping of modes
t_length = 600; % 设定计算时间长度
para_comp = six_TMD_mu_002_layout.para_comp;
TMDs_mass = para_comp(:,2);
TMDs_frequency = para_comp(:,4);
TMDs_damping_ratio = para_comp(:,6);
TMDs_location = para_comp(:,8);
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
    [result] = test_a_0_robust(number_of_modes_to_control, number_of_modes_to_consider, number_of_tmds, modal_damping_ratios, t_length, TMDs_mass, TMDs_frequency, TMDs_damping_ratio, TMDs_location,freqratio);
    S2_data(k1)=result.dis_all_modes_sum;
    dis_all_modes=result.dis_all_modes;
%     S2_data_all(:,k1)=dis_all_modes(1:6,:);
end

dis_TMD=result.dis_TMD;

t=result.t;
seqs_allmodes_max_point_all=result.seqs_allmodes_max_point_all;
% 
% subplot(3,2,1)
% plot(t,seqs_allmodes_max_point_all(1,:))
% 
% subplot(3,2,2)
% plot(t,seqs_allmodes_max_point_all(2,:))
% 
% subplot(3,2,3)
% plot(t,seqs_allmodes_max_point_all(3,:))
% 
% subplot(3,2,4)
% plot(t,seqs_allmodes_max_point_all(4,:))
% 
% subplot(3,2,5)
% plot(t,seqs_allmodes_max_point_all(5,:))
% 
% subplot(3,2,6)
% plot(t,seqs_allmodes_max_point_all(6,:))
% 
% subplot(2,2,1)
% plot(t,seqs_allmodes_max_point_all(1,:))
% 
% subplot(2,2,2)
% plot(t,dis_TMD(:,1))
% 
% subplot(2,2,3)
% plot(t,dis_TMD(:,2))
% 
% subplot(2,2,4)
% plot(t,dis_TMD(:,3))
% 
% h=t(2)-t(1);
% fs=1/h;
% time_history_cut=seqs_allmodes_max_point_all(round(end / 8 * 6, 0):end);
% 
% [psd_avg, f, psd_plot] = fft_transfer(fs,time_history_cut');
% figure
% plot(f,psd_plot)
% 
% MM = result.MM;
% CC = result.CC;
% KK = result.KK;
% 
% CC_aero=CC;
% rho=result.rho;
% U=result.U;
% a1=result.a1;
% D=result.D;
% mode_integral_2=result.mode_integral_2;
% mode_number = 1;
% CC_aero(mode_number,mode_number)=CC_aero(mode_number,mode_number)-rho*U*D*a1*mode_integral_2;
% 
% Mode_aero = Complex_Eigenvalue_Analysis(MM,CC_aero,KK)
% Mode = Complex_Eigenvalue_Analysis(MM,CC,KK)




%%

damping_ratio_ratios = 33;

%% Scheme 2
six_TMD_mu_002_layout = importdata('6TMD_mu_002_layout.mat');


modeinfo = load('modeinfo_all.mat'); % 读取未安装TMD时的模态信息. Read the mode information. 
nodegap=modeinfo.nodegap;
mode_re=modeinfo.mode_re;


nTMD = number_of_tmds;
nModes = number_of_modes_to_consider;
modal_damping_ratios = ones(1, number_of_modes_to_consider) * 0.003; % 模态阻尼 The damping of modes
% t_length = 500; % 设定计算时间长度
para_comp = six_TMD_mu_002_layout.para_comp;
TMDs_mass = para_comp(:,2);
TMDs_frequency = para_comp(:,4);
TMDs_damping_ratio = para_comp(:,6);
TMDs_location = para_comp(:,8);
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
    [result] = test_a_0_robust(number_of_modes_to_control, number_of_modes_to_consider, number_of_tmds, modal_damping_ratios, t_length, TMDs_mass, TMDs_frequency, TMDs_damping_ratio, TMDs_location,freqratio);
    S2_data(k1)=result.dis_all_modes_sum;
    dis_all_modes=result.dis_all_modes;
%     S2_data_all(:,k1)=dis_all_modes(1:6,:);
end

dis_TMD=result.dis_TMD;

t2=result.t;
seqs_allmodes_max_point_all2=result.seqs_allmodes_max_point_all;

% subplot(3,2,1)
% plot(t,seqs_allmodes_max_point_all(1,:))
% 
% subplot(3,2,2)
% plot(t,seqs_allmodes_max_point_all(2,:))
% 
% subplot(3,2,3)
% plot(t,seqs_allmodes_max_point_all(3,:))
% 
% subplot(3,2,4)
% plot(t,seqs_allmodes_max_point_all(4,:))
% 
% subplot(3,2,5)
% plot(t,seqs_allmodes_max_point_all(5,:))
% 
% subplot(3,2,6)
% plot(t,seqs_allmodes_max_point_all(6,:))
% 
% subplot(2,2,1)
% plot(t,seqs_allmodes_max_point_all(1,:))
% 
% subplot(2,2,2)
% plot(t,dis_TMD(:,1))
% 
% subplot(2,2,3)
% plot(t,dis_TMD(:,2))
% 
% subplot(2,2,4)
% plot(t,dis_TMD(:,3))
% 
% h=t(2)-t(1);
% fs=1/h;
% time_history_cut=seqs_allmodes_max_point_all(round(end / 8 * 6, 0):end);
% 
% [psd_avg, f, psd_plot] = fft_transfer(fs,time_history_cut');
% figure
% plot(f,psd_plot)
% 
% MM = result.MM;
% CC = result.CC;
% KK = result.KK;
% 
% CC_aero=CC;
% rho=result.rho;
% U=result.U;
% a1=result.a1;
% D=result.D;
% mode_integral_2=result.mode_integral_2;
% mode_number = 1;
% CC_aero(mode_number,mode_number)=CC_aero(mode_number,mode_number)-rho*U*D*a1*mode_integral_2;
% 
% Mode_aero = Complex_Eigenvalue_Analysis(MM,CC_aero,KK)
% Mode = Complex_Eigenvalue_Analysis(MM,CC,KK)
% 



%%
% figure
% subplot(3,2,1)
% plot(t,seqs_allmodes_max_point_all(1,:))
% hold on
% plot(t2,seqs_allmodes_max_point_all2(1,:))
% 
% subplot(3,2,2)
% plot(t,seqs_allmodes_max_point_all(2,:))
% hold on
% plot(t2,seqs_allmodes_max_point_all2(2,:))
% 
% subplot(3,2,3)
% plot(t,seqs_allmodes_max_point_all(3,:))
% hold on
% plot(t2,seqs_allmodes_max_point_all2(3,:))
% 
% subplot(3,2,4)
% plot(t,seqs_allmodes_max_point_all(4,:))
% hold on
% plot(t2,seqs_allmodes_max_point_all2(4,:))
% 
% subplot(3,2,5)
% plot(t,seqs_allmodes_max_point_all(5,:))
% hold on
% plot(t2,seqs_allmodes_max_point_all2(5,:))
% 
% subplot(3,2,6)
% plot(t,seqs_allmodes_max_point_all(6,:))
% hold on
% plot(t2,seqs_allmodes_max_point_all2(6,:))
% 
% 
% %%
% h=t(2)-t(1);
% fs=1/h;
% time_history_cut=seqs_allmodes_max_point_all(1,round(end / 10 * 8, 0):end);
% time_history_cut2=seqs_allmodes_max_point_all(1,round(end / 10 * 1, 0):round(end / 10 * 2, 0));
% [psd_avg, f, psd_plot] = fft_transfer(fs,time_history_cut');
% [psd_avg2, f2, psd_plot2] = fft_transfer(fs,time_history_cut2');
% figure
% plot(f,psd_plot)
% hold on 
% plot(f2,psd_plot2)


%% 
figure
plot(t,seqs_allmodes_max_point_all)
hold on
plot(t2,seqs_allmodes_max_point_all2)