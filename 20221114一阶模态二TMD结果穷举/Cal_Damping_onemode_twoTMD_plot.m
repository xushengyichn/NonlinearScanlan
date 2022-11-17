%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%Date: 2022-10-14 11:41:51
%LastEditors: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%LastEditTime: 2022-11-15 17:06:08
%FilePath: \NonlinearScanlan\20221114一阶模态二TMD结果穷举\Cal_Damping_onemode_twoTMD.m
%Description: 直接计算安装TMD后的阻尼比
%
%Copyright (c) 2022 by xushengyichn 54436848+xushengyichn@users.noreply.github.com, All Rights Reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc;clear;close all;
addpath("../函数/")


numberofTMD = 2; % 所需要计算的TMD的数量.

savedata = 0;

nTMD = numberofTMD;
% nTMD = 2; %TMD数量
modeinfo = load('modeinfo.mat');
% 设计第一个TMD的参数
% TMD1按照规范设计
fs=modeinfo.Freq(1);
mode1 = modeinfo.eig_vec(:, 1);
maxphi1=max(mode1);
phi1=maxphi1;
phi2_all=(0.001:0.001:1.5)*maxphi1;
mu = 0.02;
mTMD1 = mu / (max(mode1)^2);
fTMD1 = 1/(1+mu)*fs;
zetaTMD1 = sqrt(3*mu/8/(1+mu));

xTMD2_all=0:0.5:660;
fTMD2_all=fTMD1*0.7:0.01:fTMD1*1.4;
u_temps=0:1:0.05/maxphi1;

[U_temps,FTMD2_all]=ndgrid(u_temps,fTMD2_all);
variables = [U_temps(:),FTMD2_all(:)];
% variables = fTMD2_all';
onemode_twotmd_dis=zeros(size(variables,1),4);%第四列为是否收敛标志
numIterations=size(variables,1);


%% 只安装一个TMD的情况(模拟结构频率失调的情况)

for k2 = 1:length(u_temps)
    k1=1;
    u_temp=u_temps(k2);
    nodeondeck = importdata('nodeondeck.txt');
    KMmapping = importmappingmatrix('KMatrix.mapping');
    KK_eq = modeinfo.KK_eq(1:1, 1:1);
    KK_eqs=KK_eq*0.7:0.1:KK_eq*1.3;
    KK_eqs=KK_eq;
    eig_val=modeinfo.eig_val(1:1, 1:1);
    eig_vals = eig_val*0.7:0.1:eig_val*1.3;
    eig_vals =eig_val;
    
    mass_six_span = 10007779.7;
    mTMD = [mTMD1];
    Ftmd = [fTMD1];
    % Ftmd = [variables(k1,1) variables(k1,1)];
    zetaTMD = [zetaTMD1 ];
    % xTMD = [384.5 variables(k1,1)];
    % phiTMD = [phi1;variables(k1,1)];
    phiTMD = [phi1];
    omegaTMD = 2 * pi * Ftmd;
    
    mode_number = 1;
    ifcalmode = 3;
    h = 0.01;
    t_length =150;
    
    
    t = 0:h:t_length; % Time
    
    
    
    
        calmodes = mode_number; %考虑模态数 Consider the number of modes
        nModes = calmodes;
    
        MM_eq = modeinfo.MM_eq(1:calmodes, 1:calmodes);
        KK_eq = KK_eqs(k1);
        eig_val = eig_vals(k1);
        eig_vec = modeinfo.eig_vec(:, 1:calmodes);
    
    
          modeTMD=phiTMD;
    
        
        [result] = CalDamping_Polynomial_withTMD_multidegree_usephi_plot(1,mTMD,zetaTMD,omegaTMD,modeTMD,mode_number,ifcalmode,MM_eq,KK_eq,calmodes,eig_val,eig_vec,u_temp);
        sys_info_temp=result.Mode_sys.Mode;
        fre_info(:,k1)=sys_info_temp.Frequency;
        zeta_info(:,k1)=sys_info_temp.("Damping ratio");
    
        aero_info_temp=result.Mode_add_ADF.Mode;
        aero_fre_info(:,k1)=aero_info_temp.Frequency;
        aero_zeta_info(:,k2)=aero_info_temp.("Damping ratio");
    
        aero_complex_info(:,k1)=result.Mode_add_ADF.eig_val;
    
    
end

figure
hold on
omegas=sqrt(KK_eqs);
fs=omegas/2/pi;
plot(u_temps,aero_zeta_info(1,:))
plot(u_temps,aero_zeta_info(2,:))
% plot(fs,zeta_info)
% plot(fs,aero_zeta_info)

% scatter(fs,zeta_info)
% scatter(fs,aero_zeta_info)
singleplotdata.fs=fs;
singleplotdata.aero_zeta_info=aero_zeta_info;
% close all
% figure
% hold on
% frequency_ratio=variables;
% % frequency_ratio=variables/0.8175;
% % plot(frequency_ratio,fre_info(1,:))
% % plot(frequency_ratio,fre_info(2,:))
% % plot(frequency_ratio,fre_info(3,:))
% plot(frequency_ratio,aero_fre_info(1,:))
% plot(frequency_ratio,aero_fre_info(2,:))
% % plot(frequency_ratio,aero_fre_info(3,:))


% figure
% hold on
% plot(frequency_ratio,zeta_info(1,:))
% plot(frequency_ratio,zeta_info(2,:))
% % plot(frequency_ratio,zeta_info(3,:))
% % plot(frequency_ratio,aero_zeta_info(1,:))
% % plot(frequency_ratio,aero_zeta_info(2,:))
% % plot(frequency_ratio,aero_zeta_info(3,:))


clear fre_info zeta_info aero_fre_info aero_zeta_info aero_complex_info
%% 安装两个TMD的情况
parfor k1 = 1:size(variables,1)

mass_six_span = 10007779.7;
mTMD = [mTMD1 mTMD1];
Ftmd = [fTMD1 variables(k1,2)];
% Ftmd = [variables(k1,1) variables(k1,1)];
zetaTMD = [zetaTMD1 zetaTMD1];
% xTMD = [384.5 variables(k1,1)];
% phiTMD = [phi1;variables(k1,1)];
phiTMD = [phi1;phi1];
omegaTMD = 2 * pi * Ftmd;
u_temp=variables(k1,1);
mode_number = 1;
ifcalmode = 3;
h = 0.01;
t_length =150;


nodeondeck = importdata('nodeondeck.txt');
KMmapping = importmappingmatrix('KMatrix.mapping');
t = 0:h:t_length; % Time




    calmodes = mode_number; %考虑模态数 Consider the number of modes
    nModes = calmodes;

    MM_eq = modeinfo.MM_eq(1:calmodes, 1:calmodes);
    KK_eq = modeinfo.KK_eq(1:calmodes, 1:calmodes);
    eig_val = modeinfo.eig_val(1:calmodes, 1:calmodes);
    eig_vec = modeinfo.eig_vec(:, 1:calmodes);


      modeTMD=phiTMD;


    [result] = CalDamping_Polynomial_withTMD_multidegree_usephi_plot(nTMD,mTMD,zetaTMD,omegaTMD,modeTMD,mode_number,ifcalmode,MM_eq,KK_eq,calmodes,eig_val,eig_vec,u_temp);
    sys_info_temp=result.Mode_sys.Mode;
    fre_info(:,k1)=sys_info_temp.Frequency;
    zeta_info(:,k1)=sys_info_temp.("Damping ratio");

    aero_info_temp=result.Mode_add_ADF.Mode;
    aero_fre_info(:,k1)=aero_info_temp.Frequency;
    aero_zeta_info(:,k1)=aero_info_temp.("Damping ratio");

    aero_complex_info(:,k1)=result.Mode_add_ADF.eig_val;



end
aero_zeta_info_min=min(aero_zeta_info);
bridge_damping_grid=griddata(variables(:,1),variables(:,2),aero_zeta_info_min,U_temps,FTMD2_all);


save cal_damping_onemode_twoTMD_plot.mat 

%% 绘图
u_temps=u_temps;
singleplotdata_aero_zeta_info=singleplotdata.aero_zeta_info;
singleplotdata_min=min(singleplotdata_aero_zeta_info);
fsingle=1;
[U_temps2_single,FTMD2_all2]=ndgrid(u_temps,fTMD2_all);
variables2 = [U_temps2_single(:),FTMD2_all2(:)];
singleplotdata_min_rep=repmat(singleplotdata_min,1,58);
bridge_damping_grid_single=griddata(variables2(:,1),variables2(:,2),singleplotdata_min_rep,U_temps2_single,FTMD2_all2);
figure
surf(U_temps*maxphi1,FTMD2_all,bridge_damping_grid)
hold on
surf(U_temps2_single*maxphi1,FTMD2_all,bridge_damping_grid_single)
save alldata.mat
save damping_plot U_temps maxphi1 FTMD2_all bridge_damping_grid U_temps2_single bridge_damping_grid_single

%% 所需函数

function result = P_eq(mode, temp_vec, Matrix)
    vec = temp_vec(:, mode);
    result = vec' * Matrix * vec;
end
