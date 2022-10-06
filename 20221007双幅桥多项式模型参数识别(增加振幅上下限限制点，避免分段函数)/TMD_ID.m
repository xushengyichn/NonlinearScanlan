%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: Shengyi Xu xushengyichn@outlook.com
%Date: 2022-07-08 13:16:01
%LastEditors: Shengyi Xu xushengyichn@outlook.com
%LastEditTime: 2022-07-10 15:20:23
%FilePath: \twindeck_ID\TMD_ID.m
%Description: 该代码用来识别TMD的阻尼和频率参数
%
%Copyright (c) 2022 by Shengyi Xu xushengyichn@outlook.com, All Rights Reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear;close all;

%% 各种参数设置
recorddata = 1; % 是否需要记录数据
fitdecide =0; % 是否需要拟合振幅包络线
fname = ['SZTD-110-TMD8-1']; %分析工况对应文件名
casenum=['-2'];
fnamenew = strcat(fname,casenum);
chanum = 8; %记录通道数

%NTNU笔记本路径
addpath(genpath("C:\Users\shengyix\OneDrive\NAS云同步\Drive\0研究生\有用的代码\HHT-Tutorial-master")) %经验包络法求瞬时频率和振幅必备
addpath("C:\Users\shengyix\OneDrive\NAS云同步\Drive\0研究生\有用的代码\Matlab_PlottingTemplates") %绘图工具
%台式机路径
addpath(genpath("C:\Users\xushe\OneDrive\NAS云同步\Drive\0研究生\有用的代码\HHT-Tutorial-master"))
addpath("C:\Users\xushe\OneDrive\NAS云同步\Drive\0研究生\有用的代码\Matlab_PlottingTemplates")
addpath(genpath("D:\Matlab\tftb-0.2"))


%% 添加试验数据路径
path(1) = "D:\资料存档\实验数据\2021年4月2日深中通道110m和60m连续梁涡振试验数据\测振\20210323\TMD标定\1"; %台式机路径
path(2) =  "C:\Users\shengyix\Documents\2021年4月2日深中通道110m和60m连续梁涡振试验数据\测振\20210323\TMD标定\1";%NTNU笔记本路径

for k1 = 1:length(path)

    if exist(path(k1), 'dir') ~= 0
        addpath(path(k1))
        decide = 1;
        disp("添加试验数据路径:" + path(k1))
    end

end

if decide == 0
    error("未找到试验数据，请添加试验数据文件路径！程序终止")
end

clear decide


%% 创建响应数据分割起点与终点记录文件

%只为创建记录变量，若未找到文件会创建新文件
if exist("TMD_logfile.mat", 'file')
    disp("检测到已有数据文件，读入数据文件")
    load TMD_logfile.mat
else
      % Make N by 2 matrix of fieldname + value type
      variable_names_types = [["casename", "string"]; ...
      ["start", "double"]; ...
      ["end", "double"]; ...
      ["zeta", "double"]; ...
      ["fre", "double"]; ...
      ["f_resolution", "double"]; ...
      ];
% Make table using fieldnames & value types from above
my_table = table('Size', [0, size(variable_names_types, 1)], ...
'VariableNames', variable_names_types(:, 1), ...
'VariableTypes', variable_names_types(:, 2));

save TMD_logfile my_table;
disp("未检测到已有数据文件，创建数据文件")
end



%% 提取响应数据
for ii = 1:size(fname, 1)

    for i = 1:chanum
        refname = strcat(fname(ii, :), '#', int2str(i), '.tsp');
        fid = fopen(refname, 'r');
        fsamp = fscanf(fid, '%f,%*d,%*d,%*d,%*d,%*f,%*d,%f'); %只读取tsp文件的最后一个参数
        up_amp(i) = fsamp(2, 1);
        fs = fsamp(1, 1);
        %amp(i)为DASP中的通道传感器的标定值  amp(i)
        fclose(fid);
    end

    for j = 1:chanum
        dataname = strcat(fname(ii, :), '#', int2str(j), '.sts');
        fid = fopen(dataname, 'rb');
        [DATA(:, j), cn] = fread(fid, 'single'); %读取二进制文件，新买的DASP
        fclose(fid);
        DATA(:, j) = DATA(:, j) ./ 50; %进制转换，新买的dasp        10000mv， 10v对应量程20cm，200mm
        D(:, j) = DATA(:, j);
    end

    up_D = D(:, 1:4); %
    down_D = D(:, 5:8);

    n = length(D(:, 1));

    for i = 1:n
        AM(i, 1) = 1 / fs * i;
    end

    t = AM(:, 1); %未经过裁切的时间序列

    up_index_start=1;
    up_index_end=length(t);
    down_index_start=1;
    down_index_end=length(t);

    disp("Windward girder: Start point: " + num2str(up_index_start) + ";End point: " + num2str(up_index_end) + ";")
    disp("Leeward girder: Start point: " + num2str(down_index_start) + ";End point: " + num2str(down_index_end) + ";")
    up_D = up_D(up_index_start:up_index_end, :);
    down_D = down_D(down_index_start:down_index_end, :);
    up_t = t(up_index_start:up_index_end);
    down_t = t(down_index_start:down_index_end);
    up_t = up_t - up_t(1);
    down_t = down_t - down_t(1);

    UP(:, 1) = (up_D(:, 1) - mean(up_D(:, 1)) + up_D(:, 2) - mean(up_D(:, 2)) + up_D(:, 3) - mean(up_D(:, 3)) + up_D(:, 4) - mean(up_D(:, 4))) / 4/1000; %上游梁振动响应（m）
    DOWN(:, 1) = (down_D(:, 1) - mean(down_D(:, 1)) + down_D(:, 2) - mean(down_D(:, 2)) + down_D(:, 3) - mean(down_D(:, 3)) + down_D(:, 4) - mean(down_D(:, 4))) / 4/1000; %下游梁振动响应（m）
end
TMD_dis=down_D(:,3);%选择TMD的位移时程
TMD_t=down_t;%选择TMD的时间时程
%% 绘制响应曲线
% The standard values for colors saved in PLOT_STANDARDS() will be accessed from the variable PS
PS = PLOT_STANDARDS();
figure(1)
fig1_comps.fig = gcf;
hold on
fig1_comps.p1 = plot(TMD_t, TMD_dis);%
hold off
close all
seq1=TMD_t;

%筛选有用数据
Name = my_table.casename;
isexist = find(Name == fnamenew);


if isempty(isexist)
    disp("未检测到识别记录，请选择曲线起始点和终点。")
    decidesel = 1;
else
    showstr = "已检测到TMD识别记录，请输入是否需要重新识别（0/不需要，1/需要）。";
    decidesel=input(showstr);
end

if decidesel==0
    index_start = my_table.start(isexist);
    index_end = my_table.end(isexist);
else
    figure('Name', 'Select segment')
    plot(seq1, TMD_dis)
    disp('Please choose the segment:')
    [Xx, ~] = ginput(2);
    index_start = find(seq1 >= Xx(1), 1);
    index_end = find(seq1 >= Xx(2), 1);
end
close

TMD_t_sel = TMD_t(index_start:index_end)-TMD_t(index_start);
TMD_dis_sel = TMD_dis(index_start:index_end)-mean(TMD_dis(index_start:index_end));
plot(TMD_t_sel, TMD_dis_sel);
fs=1/mean(diff(TMD_t_sel));%差分后的采样频率（是一个采样频率的平均值，对于采样间隔存在变化的情况下是一个平均采样频率）

[~, TMD_f, TMD_plot] = fft_transfer(fs, TMD_dis_sel);
figure
plot(TMD_f, TMD_plot)
[TMD_psdmax, TMD_psdmaxseq] = max(TMD_plot);
TMD_Fre_vibration = TMD_f(TMD_psdmaxseq); %计算TMD振动频率
TMD_omega0_vibration = 2*pi*TMD_Fre_vibration; %计算TMD振动频率对应的圆频率


[TMD_yupper,~] = envelope(TMD_dis_sel, fs / 8, 'peak');
[fitresult, gof]=epsilonfit(TMD_t_sel, TMD_yupper);
Cdamping=fitresult.b;
Zeta0_TMD=Cdamping/(TMD_omega0_vibration);


disp("TMD的振动频率为："+TMD_Fre_vibration)
disp("TMD的阻尼比为："+Zeta0_TMD)

%重构TMD衰减曲线
figure;hold on;
plot(TMD_t_sel, TMD_dis_sel, 'r-')
plot(TMD_t_sel, TMD_yupper, 'r-')
plot(TMD_t_sel, fitresult.a*exp(-fitresult.b * TMD_t_sel)+fitresult.c, 'g-', 'linewidth', 1.)
limitrange=0.008*ones(length(TMD_t_sel));
plot(TMD_t_sel,limitrange);
plot(TMD_t_sel,-limitrange);
legend('Raw Data', 'Amplitude Data', 'Attenuation Curve','Precision limit')
ylabel("displacement/mm")
xlabel("t/s")


f_resolution=1/TMD_t_sel(end);%频率分辨率

%去除大振幅下拟合偏移区间
front_end=0.3;


% % 短时傅里叶变换
% N= length(TMD_dis_sel);%采样点数
% fs=1/mean(diff(TMD_t_sel));%差分后的采样频率（是一个采样频率的平均值，对于采样间隔存在变化的情况下是一个平均采样频率）
% win=hanning(255);%窗函数
% [B,t,f]=tfrstft(TMD_dis_sel,1:N,N,win);%短时傅里叶变换
% figure
% imagesc(TMD_t_sel,f(1:N/2)*fs,abs(B(1:N/2,:))); axis xy
% ylim([0 10])
 if recorddata ==1
     T.casename=fnamenew;
     T.start=index_start;
     T.end=index_end;
     T.zeta=Zeta0_TMD;
     T.fre=TMD_Fre_vibration;
     T.f_resolution=f_resolution;

     Name = my_table.casename;
     isexist = find(Name == fnamenew);
    if isempty(isexist)
        alertstr = "该TMD尚未含有数据，是否添加该工况数据？（0/不添加；1/添加）";
        decidedata = input(alertstr);
        if decidedata == 1
            my_table_new=struct2table(T);
            my_table = [my_table; my_table_new];
            save TMD_logfile my_table;
            disp(T.casename + "工况已经添加")
        else
            if decidedata == 0
                disp(T.casename + "工况未添加")
            else
                disp("输入错误，工况未添加")
            end
        end
    else
        alertstr = "该数据已经存在，是否需要更新？（0/不更新；1/更新）";
        decidedata = input(alertstr);

        if decidedata == 1
            my_table_new=struct2table(T);
            my_table(isexist, :) = my_table_new;
            save TMD_logfile my_table;
            disp(T.casename + "工况已经更新")
        else

            if decidedata == 0
                disp(T.casename + "工况未更新")
            else
                disp("输入错误，工况未更新")
            end

        end

    end
 end

close all



 function [fitresult, gof] = epsilonfit(TMD_t_sel, TMD_yupper)
%CREATEFIT1(TMD_T_SEL,TMD_YUPPER)
%  Create a fit.
%
%  Data for 'epsilonfit' fit:
%      X Input : TMD_t_sel
%      Y Output: TMD_yupper
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 10-Jul-2022 15:03:20


%% Fit: 'epsilonfit'.
[xData, yData] = prepareCurveData( TMD_t_sel, TMD_yupper );

% Set up fittype and options.
ft = fittype( 'a*exp(-b*x)+c', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.240250585015601 0.180502599740628 0.844048426282812];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'epsilonfit' );
h = plot( fitresult, xData, yData );
legend( h, 'TMD_yupper vs. TMD_t_sel', 'epsilonfit', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'TMD_t_sel', 'Interpreter', 'none' );
ylabel( 'TMD_yupper', 'Interpreter', 'none' );
grid on
end



