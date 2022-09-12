%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: Shengyi Xu xushengyichn@outlook.com
%Date: 2022-06-27 16:21:36
%LastEditors: Shengyi Xu xushengyichn@outlook.com
%LastEditTime: 2022-09-12 23:25:44
%FilePath: \NonlinearScanlan\test_Polynomial_withTMD_singledegree.m
%Description: 本代码是用于求解多项式模型下，测试节段模型安装tmd后产生多阶模态，最后计算得到响应的频率问题。
%
%Copyright (c) 2022 by Shengyi Xu xushengyichn@outlook.com, All Rights Reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clc
clear
% This is a function
% function out = Polynomial_withTMD_singledegree(Fre, Mass, Zeta0, rho, D, U, Y1k, epsilonk, Y2k, ClK, t, P, u0, udot0, mtmd, ctmd, ktmd)

% Nonlinear Newmark's Direct Integration Method with polynomial model
% (n = number of time steps)
% (ndof = number degrees of freedom)

% INPUT
% Fre      = Frequency of the system         =>[1]
% Mass     = Mass of the system              =>[1]
% Zeta0    = Damping ratio of the system     =>[1]
% rho      = Air density                     =>[1]
% D        = Reference length                =>[1]
% U        = Wind speed at a certain reduced frequency       =>[1]
% Y1k      = Y1 at a certain reduced frequency       =>[1]
% epsilonk = Epsilon at a certain reduced frequency       =>[1]
% Y2k      = Y2 at a certain reduced frequency       =>[1]
% Clk      = Cl at a certain reduced frequency       =>[1]
% t        = Time vector         =>[1,n]
% P        = load vs. time       =>[1,n]
% u0       = Initial displacements =>[1]
% udot0    = Initial velocity =>[1]
% gam      = gamma (constant)
% beta     = beta  (constant)
% mtmd     = mass of tmd
% ctmd     = damping of tmd
% ktmd     = stiffness of tmd

%--------------------------------------------------------------------------
% beta = 0,     gamma = 1/2 => explicit central difference method
% beta = 1/4,   gamma = 1/2 => undamped trapezoidal rule (implicit)

%--------------------------------------------------------------------------
gamma = 1/2; % Factor in the Newmark algorithm
beta = 1/4; % Factor in the Newmark algorithm
%--------------------------------------------------------------------------
% 设置矩阵大小
matrixsize = 5; %calculate one degree of freedom
%--------------------------------------------------------------------------
% 导入试验数据
load('TMD_logfile.mat');
my_table_tmd=my_table;
load('SZTD110_logfile.mat');

%% 提取工况设置
fname = ['SZTD-110-case2-22.3-fasan-2401']; %记录文件名

chanum = 8; %记录通道数
filename = strsplit(fname, '-');
casenumber = cell2mat(filename(3));
spacing = str2double(cell2mat(filename(4)));
type = cell2mat(filename(5));
Rspeed = cell2mat(filename(6));

casenumber = str2double(casenumber(5:end));
Rspeed = Rspeed(1:end - 1);
Rspeed = str2double(Rspeed);

if strcmp(type, 'fasan')
    type = 1;
else

    if strcmp(type, 'shuaijian')
        type = 2;
    else
        error("实验数据文件名可能有误，请确认信号状态为衰减或发散，程序终止")
    end

end

% 判断工况风攻角
if casenumber <= 17 || and(casenumber >= 32, casenumber <= 35) || and(casenumber >= 40, casenumber <= 16) || casenumber == 55
    AOA = 3;
else

    if and(casenumber >= 18, casenumber <= 25) || and(casenumber >= 36, casenumber <= 37) || and(casenumber >= 47, casenumber <= 50) || casenumber == 54
        AOA = 0;
    else

        if and(casenumber >= 26, casenumber <= 31) || and(casenumber >= 38, casenumber <= 39) || and(casenumber >= 51, casenumber <= 53) || and(casenumber >= 56, casenumber <= 57)
            AOA = -3;
        end

    end

end

%导入试验数据
Name = my_table.casename;
isexist = find(Name == fname);
up_a1=my_table.up_parameter_a1(isexist);
up_a2=my_table.up_parameter_a2(isexist);
up_a3=my_table.up_parameter_a3(isexist);
up_a4=my_table.up_parameter_a4(isexist);
up_a5=my_table.up_parameter_a5(isexist);

up_a = [up_a1 up_a2 up_a3 up_a4 up_a5];

up_H4=my_table.up_parameter_H4(isexist);% 气动刚度
up_upperlimit= my_table.up_upperlimit(isexist); %除以特征长度D的无量纲振幅
up_lowerlimit= my_table.up_lowerlimit(isexist); %除以特征长度D的无量纲振幅
% up_upperlimit=10000;
% up_lowerlimit=-10000;
up_Fren_vibration_withwind=my_table.up_Fren_vibration_withwind(isexist); 


% 结构参数
% Structural parameters
D = 0.667; % deck depth
m = 80; % mass of the segment model
Mass = m;
F0 = my_table.up_Fre_vibration(isexist); % Frequency without wind
Fre= F0;
rho = 1.225; % density of the air
U = my_table.Windspeed(isexist); % density of the air
% U = 5;
Zeta0 =  my_table.up_dltx_zeta0(isexist); % damping ratio without wind

h = 1/256; % 时间步长 % Step size of the algorithm
up_t=0:h:60; % Time vector
up_told = up_t;
up_tt = up_told;
up_P = zeros(1, size(up_told, 1));

disp("节段模型质量："+num2str(Mass));
disp("节段模型阻尼系数："+num2str(Zeta0));
disp("节段模型频率："+num2str(Fre));
disp("节段模型风速："+num2str(U));
disp("节段模型气动阻尼参数"+num2str(up_a));
disp("节段模型气动刚度参数"+num2str(up_H4));


%% TMD 参数
sel=[12 14 15 18];
% sel=[4 7 8 10];
% sel=[11 11 11 11];
% sel=[11];
zetatmd = my_table_tmd.zeta(sel);
fretmd = my_table_tmd.fre(sel)+0.16;
mtmd = ones(length(sel),1)*0.25;
disp("TMD质量："+num2str(mtmd));
disp("TMD阻尼系数："+num2str(zetatmd));
disp("TMD频率："+num2str(fretmd));
%% Calculate the response


up_u0 = [-0.086; -1;-1;-1;-1]*10e-3;
up_udot0 = [0; 0;0;0;0];

P = zeros(5, length(up_tt));
nModes = 1;
matrixsize=5;

% 
% up_u0 = [-0.086; -1]*10e-3;
% up_u0 = [-0.14534; 1]*10e-3;
% up_udot0 = [0; 0];
% 
% P = zeros(2, length(up_tt));
% nModes = 1;
% matrixsize=2;

[MM,CC,KK]=CreateMatrixwithTMD(nModes,Mass,Zeta0,Fre,mtmd,zetatmd,fretmd);

a1_lower=up_a1+4/3*up_a2/pi*up_lowerlimit+up_a3/4*up_lowerlimit^2+8/15*up_a4/pi*up_lowerlimit^3+up_a5/8*up_lowerlimit^4;
disp("lowerlimit="+a1_lower)
% a1_lower=a1_lower*2
[CC1,KK1]=AddAerodynamicDampingandStiffness(CC,KK,rho,U,D,a1_lower,up_H4);

[V,DD]=eigs(KK1,MM);
Result0=sort(diag(sqrt(DD)/2/pi));
Result1= Complex_Eigenvalue_Analysis(MM,CC,KK1);%不考虑气动阻尼的复数特征值分析
Result2=Complex_Eigenvalue_Analysis(MM,CC1,KK1);%考虑气动阻尼的复特征值分析
disp(Result0)
% up_a(1)=a1_lower;
% up_a(2)=0;
% up_a(3)=0;
% up_a(4)=0;
% up_a(5)=0;
out = test_polynomial_NB_withTMDs_addstiff_withlimit(Fre, Mass, Zeta0, rho, D, U, up_a,up_H4, up_t,h, P, up_u0, up_udot0,up_upperlimit,up_lowerlimit,up_Fren_vibration_withwind,nModes,mtmd,fretmd,zetatmd);

figure
subplot(1, 2, 1)
plot(out(:, 1), out(:, 2))
ylim([-0.01 0.01])
title("main structure")
xlabel('time(s)');
ylabel('displacement(m)');
subplot(1, 2, 2)
plot(out(:, 1), out(:, 3))
title("TMD")
ylim([-0.01 0.01])

fs=1/(up_t(2)-up_t(1));
caldata=out(round(end/2,0):end,2);
[psd_avg, f2, psd_plot2] = fft_transfer(256,caldata);
figure
plot(f2, psd_plot2)
xlabel("Frequency")
title("计算响应的频谱")
disp("主结构响应均方根值"+num2str(std(out(:, 2))))
[a,b]=max(psd_plot2);
disp("振动频率为："+num2str(f2(b)))
disp(Result1)
disp(Result2)
% % figure 
% % plot(out(:, 1), out(:, 2))
% % ylim([-0.01 0.01]/5)
% % title("main structure")
% % hold on
% % load test.mat
% % plot(up_t,UP)
% % legend("calculated","windtunnel test")
% fs=1/(up_t(2)-up_t(1));
% % [psd_avg, f1, psd_plot1] = fft_transfer(fs,UP);
% [psd_avg, f2, psd_plot2] = fft_transfer(256,out(:, 2));
% % figure
% % plot(f1, psd_plot1)
% % hold on 
% figure
% plot(f2, psd_plot2)
% xlabel("Frequency")
% title("计算响应的频谱")
% disp("主结构响应最大值"+num2str(max(out(:, 2))))
% % [a,b]=max(psd_plot1);
% % f1(b)
% 
% [a,b]=max(psd_plot2);
% disp("振动频率为："+num2str(f2(b)))

% 变化tmd阻尼
exe=0;
if exe == 1
zetatmds=0.001:0.001:0.2;
numIterations=length(zetatmds);
ppm = ParforProgressbar(numIterations);
freqvibration=zeros(length(zetatmds),1); %节段模型振动频率
freqTMD=zeros(length(zetatmds),1); %TMD振动频率
Freqsystem=zeros(length(zetatmds),2); 
parfor i=1:length(zetatmds)
    [out,Freq] = test_polynomial_NB_withTMDs_addstiff_withlimit(Fre, Mass, Zeta0, rho, D, U, up_a,up_H4, up_t,h, P, up_u0, up_udot0,up_upperlimit,up_lowerlimit,up_Fren_vibration_withwind,nModes,mtmd,fretmd,zetatmds(i));
    outclip=out(round(end/2,0):end,2);
    outcliptmd=out(round(end/2,0):end,3);
    [psd_avg, f2, psd_plot2] = fft_transfer(256,outclip);
    [psd_avg, ftmd, psd_plot3] = fft_transfer(256,outcliptmd);
    [a,b]=max(psd_plot2);
    [c,d]=max(psd_plot3);
    % disp("振动频率为："+num2str(f2(b)))
    freqvibration(i)=f2(b);
    freqTMD(i)=ftmd(d);
    Freqsystem(i,:)=Freq; 
    % progressbar(i/allcasenumber)
    ppm.increment();
end
delete(ppm);
figure
plot(zetatmds,freqvibration)
hold on 
plot(zetatmds,freqTMD)
plot(zetatmds,Freqsystem(:,1))
plot(zetatmds,Freqsystem(:,2))
title("Frequency varies with zeta")
legend("vibration frequency of sectional model","vibration frequency of TMD","mode frequency 1","mode frequency 2")
xlabel("zeta")
ylabel("frequency")
end


% 变化tmd频率（通过变化质量实现）
exe =0 ;
if exe ==1
mtmds=0.001:0.001:5;
numIterations=length(mtmds);
ppm = ParforProgressbar(numIterations);
freqvibration=zeros(numIterations,1); %节段模型振动频率
freqTMD=zeros(numIterations,1); %TMD振动频率
Freqsystem=zeros(numIterations,2); 
parfor i=1:numIterations
    mtmd=mtmds(i);
    [out,Freq] = test_polynomial_NB_withTMDs_addstiff_withlimit(Fre, Mass, Zeta0, rho, D, U, up_a,up_H4, up_t,h, P, up_u0, up_udot0,up_upperlimit,up_lowerlimit,up_Fren_vibration_withwind,nModes,mtmd,fretmd,zetatmd);
    outclip=out(round(end/2,0):end,2);
    outcliptmd=out(round(end/2,0):end,3);
    [psd_avg, f2, psd_plot2] = fft_transfer(256,outclip);
    [psd_avg, ftmd, psd_plot3] = fft_transfer(256,outcliptmd);
    [a,b]=max(psd_plot2);
    [c,d]=max(psd_plot3);
    % disp("振动频率为："+num2str(f2(b)))
    freqvibration(i)=f2(b);
    freqTMD(i)=ftmd(d);
    Freqsystem(i,:)=Freq; 
    % progressbar(i/allcasenumber)
    ppm.increment();
end
delete(ppm);
figure
plot(mtmds,freqvibration)
hold on 
plot(mtmds,freqTMD)
plot(mtmds,Freqsystem(:,1))
plot(mtmds,Freqsystem(:,2))
title("Frequency varies with mass of tmd")
legend("vibration frequency of sectional model","vibration frequency of TMD","mode frequency 1","mode frequency 2")
xlabel("mass of tmd")
ylabel("frequency")

end

% 变化tmd频率（通过变化刚度实现）
exe =0;
if exe ==1
    fretmds = 3:0.01:7;
    numIterations=length(fretmds);
    ppm = ParforProgressbar(numIterations);
    freqvibration=zeros(numIterations,1); %节段模型振动频率
    freqTMD=zeros(numIterations,1); %TMD振动频率
    Freqsystem=zeros(numIterations,2); 
    parfor i=1:numIterations
        fretmd=fretmds(i);
        [out,Freq] = test_polynomial_NB_withTMDs_addstiff_withlimit(Fre, Mass, Zeta0, rho, D, U, up_a,up_H4, up_t,h, P, up_u0, up_udot0,up_upperlimit,up_lowerlimit,up_Fren_vibration_withwind,nModes,mtmd,fretmd,zetatmd);
        outclip=out(round(end/2,0):end,2);
        outcliptmd=out(round(end/2,0):end,3);
        [psd_avg, f2, psd_plot2] = fft_transfer(256,outclip);
        [psd_avg, ftmd, psd_plot3] = fft_transfer(256,outcliptmd);
        [a,b]=max(psd_plot2);
        [c,d]=max(psd_plot3);
        % disp("振动频率为："+num2str(f2(b)))
        freqvibration(i)=f2(b);
        freqTMD(i)=ftmd(d);
        Freqsystem(i,:)=Freq; 
        % progressbar(i/allcasenumber)
        ppm.increment();
    end
    delete(ppm);
    figure
    plot(fretmds,freqvibration)
    hold on 
    plot(fretmds,freqTMD)
    plot(fretmds,Freqsystem(:,1))
    plot(fretmds,Freqsystem(:,2))
    title("Frequency varies with frequency of tmd")
    legend("vibration frequency of sectional model","vibration frequency of TMD","mode frequency 1","mode frequency 2")
    xlabel("frequency of tmd")
    ylabel("frequency")
end

% figure 
% plot(out(:, 1), out(:, 2))
% hold on
% plot(out1(:, 1), out1(:, 2))

% 
% % animation
% % figure
% displacement=out(:,2:6);
% % displacement(:,1)=displacement(:,1)/max(abs(displacement(:,1)));
% % displacement(:,2)=displacement(:,2)/max(abs(displacement(:,2)));
% % displacement(:,3)=displacement(:,3)/max(abs(displacement(:,3)));
% % displacement(:,4)=displacement(:,4)/max(abs(displacement(:,4)));
% % displacement(:,5)=displacement(:,5)/max(abs(displacement(:,5)));
% line=[1 2 3 4 5];
% for i = 1:1500
%     scatter(line,displacement(i+10000,:))
%     ylim([-max(abs(displacement(:,2))) max(abs(displacement(:,2)))])
%     M(:, i) = getframe;
% end
% % movie(M)
% 
% v = VideoWriter('newfile1.avi');
% open(v)
% writeVideo(v,M)
% close(v)