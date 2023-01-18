%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: Shengyi Xu xushengyichn@outlook.com
%Date: 2022-06-16 11:37:30
%LastEditors: Shengyi Xu xushengyichn@outlook.com
%LastEditTime: 2022-09-30 23:17:36
%FilePath: \twindeck_ID\twindeck_ID.m
%Description: 提取双幅桥振动响应，并使用多项式模型进行参数识别

%Copyright (c) 2022 by Shengyi Xu xushengyichn@outlook.com, All Rights Reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;
addpath('..\函数\')
%% 各种参数设置
recorddata = 0; % 是否需要记录数据
fitdecide =0; % 是否需要拟合振幅包络线
fname = ['SZTD-110-case2-22.3-fasan-2501']; %分析工况对应文件名
fname_up_dltx= ['SZTD-110-case2-22.3-dltx2'];%上游动力特性文件名
fname_down_dltx= ['SZTD-110-case2-22.3-dltx4'];%下游动力特性文件名
L = 3.6; % 节段模型长度
m_up=80/L;%上游梁单位长度质量
m_down=80/L;%下游梁单位长度质量
% m_up=80;%上游梁长度质量（虽然使用整体节段模型质量计算振幅与阻尼比关系不同，但是对气动力模型参数有影响）
% m_down=80;%下游梁长度质量
%NTNU笔记本路径
addpath(genpath("C:\Users\shengyix\OneDrive\NAS云同步\Drive\0研究生\有用的代码\HHT-Tutorial-master")) %经验包络法求瞬时频率和振幅必备
addpath("C:\Users\shengyix\OneDrive\NAS云同步\Drive\0研究生\有用的代码\Matlab_PlottingTemplates") %绘图工具
%台式机路径
addpath(genpath("C:\Users\xushe\OneDrive\NAS云同步\Drive\0研究生\有用的代码\HHT-Tutorial-master"))
addpath("C:\Users\xushe\OneDrive\NAS云同步\Drive\0研究生\有用的代码\Matlab_PlottingTemplates")
%% 添加试验数据路径
path(1) = "C:\Users\shengyix\Documents\2021年4月2日深中通道110m和60m连续梁涡振试验数据\测振\20210323\1.开槽间距6.7m\2.有栏杆\+3\2.涡振期间气动力参数识别"; %NTNU笔记本路径
path(2) = "D:\资料存档\实验数据\2021年4月2日深中通道110m和60m连续梁涡振试验数据\测振\20210323\1.开槽间距6.7m\2.有栏杆\+3\2.涡振期间气动力参数识别"; %台式机路径
path(3) = "D:\资料存档\实验数据\2021年4月2日深中通道110m和60m连续梁涡振试验数据\测振\20210323\1.开槽间距6.7m\2.有栏杆\+3\dltx"; %动力特性识别
path(4) = "C:\Users\shengyix\Documents\2021年4月2日深中通道110m和60m连续梁涡振试验数据\测振\20210323\1.开槽间距6.7m\2.有栏杆\+3\dltx"; %NTNU笔记本路径 动力特性
decide = 0;

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
if exist("SZTD110_logfile.mat", 'file')
    disp("检测到已有数据文件，读入数据文件")
    load SZTD110_logfile.mat
else
    % Make N by 2 matrix of fieldname + value type
    variable_names_types = [["casename", "string"]; ...
                            ["case", "double"]; ...
                            ["AOA", "double"]; ...
                            ["spacing", "double"]; ...
                            ["type", "double"]; ...
                            ["Rspeed", "double"]; ...
                            ["Windspeed", "double"]; ...
                            ["up_Fre_vibration", "double"]; ...
                            ["up_dltx_zeta0", "double"]; ...
                            ["up_Fren_vibration_withwind", "double"]; ...
                            ["up_ReducedFre", "double"]; ...
                            ["down_Fre_vibration", "double"]; ...
                            ["down_dltx_zeta0", "double"]; ...
                            ["down_Fren_vibration_withwind", "double"]; ...
                            ["down_ReducedFre", "double"]; ...
                            ["up_start", "double"]; ...
                            ["up_end", "double"]; ...
                            ["down_start", "double"]; ...
                            ["down_end", "double"]; ...
                            ["up_parameter_a1", "double"]; ...
                            ["up_parameter_a2", "double"]; ...
                            ["up_parameter_a3", "double"]; ...
                            ["up_parameter_a4", "double"]; ...
                            ["up_parameter_a5", "double"]; ...
                            ["up_parameter_H4", "double"]; ...
                            ["down_parameter_a1", "double"]; ...
                            ["down_parameter_a2", "double"]; ...
                            ["down_parameter_a3", "double"]; ...
                            ["down_parameter_a4", "double"]; ...
                            ["down_parameter_a5", "double"]; ...
                            ["down_parameter_H4", "double"]; ...
                            ["up_dltx_start", "double"]; ...
                            ["up_dltx_end", "double"]; ...
                            ["up_upperlimit", "double"]; ...
                            ["up_lowerlimit", "double"]; ...
                            ["down_dltx_start", "double"]; ...
                            ["down_dltx_end", "double"]; ...
                            ["down_upperlimit", "double"]; ...
                            ["down_lowerlimit", "double"]; ...
                            ];
    % Make table using fieldnames & value types from above
    my_table = table('Size', [0, size(variable_names_types, 1)], ...
    'VariableNames', variable_names_types(:, 1), ...
        'VariableTypes', variable_names_types(:, 2));

    save SZTD110_logfile my_table;
    disp("未检测到已有数据文件，创建数据文件")
end

%% 提取工况设置

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
    %筛选有用数据
    Name = my_table.casename;
    isexist = find(Name == fname);

    if isempty(isexist)
        disp("未检测到识别记录，请选择曲线起始点和终点。")
        decidesel = 1;
    else
        decidesel = 0;
        showstr = "已检测到上游梁识别记录，请输入是否需要重新识别（0/不需要，1/需要）。";
        decidesel1 = input(showstr);
        showstr = "已检测到下游梁识别记录，请输入是否需要重新识别（0/不需要，1/需要）。";
        decidesel2 = input(showstr);
    end

    seq = 1:size(up_D, 1); %上下游梁可以通用seq
    fs = 256; %激光位移计的采样频率 % Sampling frequency of the laser displacement meter

    if decidesel == 1
        figure('Name', 'Select segment')
        plot(seq, up_D(:, 1))
        disp('Please choose the segment:')
        [Xx, ~] = ginput(2);
        up_index_start = find(seq >= Xx(1), 1);
        up_index_end = find(seq >= Xx(2), 1);
        close


        disp("下游梁时程选择是否与上游梁一致？（0/一致，1/不一致）")
        decidesel_down = input("请输入：");
%         decidesel_down = 1;

        if decidesel_down == 0
            down_index_start = up_index_start;
            down_index_end = up_index_end;
        else
            figure('Name', 'Select segment')
            plot(seq, down_D(:, 1))
            disp('Please choose the segment:')
            [Xx, ~] = ginput(2);
            down_index_start = find(seq >= Xx(1), 1);
            down_index_end = find(seq >= Xx(2), 1);
            close

        end

    else

        if decidesel1 == 0
            up_index_start = my_table.up_start(isexist);
            up_index_end = my_table.up_end(isexist);
        else
            figure('Name', 'Select segment')
            plot(seq, up_D(:, 1))
            disp('Please choose the segment:')
            [Xx, ~] = ginput(2);
            up_index_start = find(seq >= Xx(1), 1);
            up_index_end = find(seq >= Xx(2), 1);
            close
        end

        if decidesel2 == 0
            down_index_start = my_table.down_start(isexist);
            down_index_end = my_table.down_end(isexist);
        else
            figure('Name', 'Select segment')
            plot(seq, down_D(:, 1))
            disp('Please choose the segment:')
            [Xx, ~] = ginput(2);
            down_index_start = find(seq >= Xx(1), 1);
            down_index_end = find(seq >= Xx(2), 1);
            close
        end

    end

    n = length(D(:, 1));

    for i = 1:n
        AM(i, 1) = 1 / fs * i;
    end

    t = AM(:, 1); %未经过裁切的时间序列

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
clear refname fid fsamp
% 提取无风状态响应数据，识别阻尼比和频率

% 读取上游梁数据
for i = 1:chanum
    refname = strcat(fname_up_dltx(ii, :), '#', int2str(i), '.tsp');
    fid = fopen(refname, 'r');
    fsamp = fscanf(fid, '%f,%*d,%*d,%*d,%*d,%*f,%*d,%f'); %只读取tsp文件的最后一个参数
    up_amp_dltx(i) = fsamp(2, 1);
    fs = fsamp(1, 1);
    %amp(i)为DASP中的通道传感器的标定值  amp(i)
    fclose(fid);
end

for j = 1:chanum
    dataname = strcat(fname_up_dltx(ii, :), '#', int2str(j), '.sts');
    fid = fopen(dataname, 'rb');
    [DATA_dltx(:, j), cn] = fread(fid, 'single'); %读取二进制文件，新买的DASP
    fclose(fid);
    DATA_dltx(:, j) = DATA_dltx(:, j) ./ 50; %进制转换，新买的dasp        10000mv， 10v对应量程20cm，200mm
    D_up_dltx(:, j) = DATA_dltx(:, j);
end
clear DATA_dltx
% 读取下游梁数据
for i = 1:chanum
    refname = strcat(fname_down_dltx(ii, :), '#', int2str(i), '.tsp');
    fid = fopen(refname, 'r');
    fsamp = fscanf(fid, '%f,%*d,%*d,%*d,%*d,%*f,%*d,%f'); %只读取tsp文件的最后一个参数
    up_amp_dltx(i) = fsamp(2, 1);
    fs = fsamp(1, 1);
    %amp(i)为DASP中的通道传感器的标定值  amp(i)
    fclose(fid);
end

for j = 1:chanum
    dataname = strcat(fname_down_dltx(ii, :), '#', int2str(j), '.sts');
    fid = fopen(dataname, 'rb');
    [DATA_dltx(:, j), cn] = fread(fid, 'single'); %读取二进制文件，新买的DASP
    fclose(fid);
    DATA_dltx(:, j) = DATA_dltx(:, j) ./ 50; %进制转换，新买的dasp        10000mv， 10v对应量程20cm，200mm
    D_down_dltx(:, j) = DATA_dltx(:, j);
end


up_D_dltx = D_up_dltx(:, 1:4); %
down_D_dltx = D_down_dltx(:, 5:8);
%筛选有用数据
Name = my_table.casename;
isexist = find(Name == fname);

if isempty(isexist)
    disp("未检测到识别记录，请选择曲线起始点和终点。")
    decidesel = 1;
else
    decidesel = 0;
    showstr = "已检测到上游梁识别记录，请输入是否需要重新识别动力特性（0/不需要，1/需要）。";
    decidesel1 = input(showstr);
    showstr = "已检测到下游梁识别记录，请输入是否需要重新识别动力特性（0/不需要，1/需要）。";
    decidesel2 = input(showstr);
end

seq1 = 1:size(up_D_dltx, 1); 
seq2 = 1:size(down_D_dltx, 1); 

if decidesel == 1
    alertstr="是否需要采用其他工况的动力特性数据结果，可以直接输入对应数据编号（不太清楚就直接输入1）,需要手动点选输入0";
    decide=input(alertstr);
    if decide~=0
        up_index_start_dltx = my_table.up_dltx_start(decide);
        up_index_end_dltx = my_table.up_dltx_end(decide);
    else
        figure('Name', 'Select segment')
        plot(seq1, up_D_dltx(:, 1))
        disp('Please choose the segment:')
        [Xx, ~] = ginput(2);
        up_index_start_dltx = find(seq1 >= Xx(1), 1);
        up_index_end_dltx = find(seq1 >= Xx(2), 1);
    end

    close

    if decide~=0
        down_index_start_dltx = my_table.down_dltx_start(decide);
        down_index_end_dltx = my_table.down_dltx_end(decide);
    else
        figure('Name', 'Select segment')
        plot(seq2, down_D_dltx(:, 1))
        disp('Please choose the segment:')
        [Xx, ~] = ginput(2);
        down_index_start_dltx = find(seq2 >= Xx(1), 1);
        down_index_end_dltx = find(seq2 >= Xx(2), 1);
    close
    end


else

    if decidesel1 == 0
        up_index_start_dltx = my_table.up_dltx_start(isexist);
        up_index_end_dltx = my_table.up_dltx_end(isexist);
    else
        figure('Name', 'Select segment')
        plot(seq1, up_D_dltx(:, 1))
        disp('Please choose the segment:')
        [Xx, ~] = ginput(2);
        up_index_start_dltx = find(seq1 >= Xx(1), 1);
        up_index_end_dltx = find(seq1 >= Xx(2), 1);
        close
    end

    if decidesel2 == 0
        down_index_start_dltx = my_table.down_dltx_start(isexist);
        down_index_end_dltx = my_table.down_dltx_end(isexist);
    else
        figure('Name', 'Select segment')
        plot(seq2, down_D_dltx(:, 1))
        disp('Please choose the segment:')
        [Xx, ~] = ginput(2);
        down_index_start_dltx = find(seq2 >= Xx(1), 1);
        down_index_end_dltx = find(seq2 >= Xx(2), 1);
        close
    end
end

n_up_dltx = length(up_D_dltx(:, 1));
n_down_dltx = length(down_D_dltx(:, 1));

for i = 1:n_up_dltx
    AM_up_dltx(i, 1) = 1 / fs * i;
end

for i = 1:n_down_dltx
    AM_down_dltx(i, 1) = 1 / fs * i;
end

t_up_dltx = AM_up_dltx(:, 1); %未经过裁切的时间序列
t_downs_dltx = AM_down_dltx(:, 1); %未经过裁切的时间序列

disp("Windward girder: Start point: " + num2str(up_index_start_dltx) + ";End point: " + num2str(up_index_end_dltx) + ";")
disp("Leeward girder: Start point: " + num2str(down_index_start_dltx) + ";End point: " + num2str(down_index_end_dltx) + ";")
up_D_dltx = up_D_dltx(up_index_start_dltx:up_index_end_dltx, :);
down_D_dltx = down_D_dltx(down_index_start_dltx:down_index_end_dltx, :);
up_t_dltx = t_up_dltx(up_index_start_dltx:up_index_end_dltx);
down_t_dltx = t_downs_dltx(down_index_start_dltx:down_index_end_dltx);
up_t_dltx = up_t_dltx - up_t_dltx(1);
down_t_dltx = down_t_dltx - down_t_dltx(1);

UP_dltx(:, 1) = (up_D_dltx(:, 1) - mean(up_D_dltx(:, 1)) + up_D_dltx(:, 2) - mean(up_D_dltx(:, 2)) + up_D_dltx(:, 3) - mean(up_D_dltx(:, 3)) + up_D_dltx(:, 4) - mean(up_D_dltx(:, 4))) / 4/1000; %上游梁振动响应（m）
DOWN_dltx(:, 1) = (down_D_dltx(:, 1) - mean(down_D_dltx(:, 1)) + down_D_dltx(:, 2) - mean(down_D_dltx(:, 2)) + down_D_dltx(:, 3) - mean(down_D_dltx(:, 3)) + down_D_dltx(:, 4) - mean(down_D_dltx(:, 4))) / 4/1000; %下游梁振动响应（m）

% 识别上游梁频率和阻尼比
[up_psd_avg_dltx, up_f_dltx, up_psd_plot_dltx] = fft_transfer(fs, UP_dltx);
[up_psdmax_dltx, up_psdmaxseq_dltx] = max(up_psd_plot_dltx);
up_Fre_vibration = up_f_dltx(up_psdmaxseq_dltx); %计算系统振动频率
up_omega0_vibration = 2*pi*up_Fre_vibration; %计算系统振动频率对应的圆频率

[up_yupper,~] = envelope(UP_dltx, fs / 8, 'peak');
up_logYupper = log(up_yupper);
p = polyfit(up_t_dltx, up_logYupper, 1);
Cdamping = -p(1);
Zeta0_up=Cdamping/(up_omega0_vibration);

% 
% figure
% title(['epsilon attenuation function']);
% plot(up_t_dltx, UP_dltx, 'r-')
% hold on
% plot(up_t_dltx, up_yupper, 'r-')
% plot(up_t_dltx, exp(p(2))*exp(p(1) * up_t_dltx), 'g-', 'linewidth', 1.)
% legend('Raw Data', 'Amplitude Data', 'Attenuation Curve')

% 识别下游梁频率和阻尼比
[down_psd_avg_dltx, down_f_dltx, down_psd_plot_dltx] = fft_transfer(fs, DOWN_dltx);
[down_psdmax_dltx, down_psdmaxseq_dltx] = max(down_psd_plot_dltx);
down_Fre_vibration = down_f_dltx(down_psdmaxseq_dltx); %计算系统振动频率
down_omega0_vibration = 2*pi*down_Fre_vibration; %计算系统振动频率对应的圆频率

[down_ydownper,~] = envelope(DOWN_dltx, fs / 8, 'peak');
down_logYdownper = log(down_ydownper);
p = polyfit(down_t_dltx, down_logYdownper, 1);
Cdamping = -p(1);
Zeta0_down=Cdamping/(down_omega0_vibration);

%% 绘制响应曲线
% The standard values for colors saved in PLOT_STANDARDS() will be accessed from the variable PS
PS = PLOT_STANDARDS();
figure(1)
fig1_comps.fig = gcf;
hold on
fig1_comps.p1 = plot(up_t, UP);
fig1_comps.p2 = plot(down_t, DOWN);
hold off

%========================================================
% ADD LABELS, TITLE, LEGEND
title('Response of the twin girders');
xlabel('Time/s');
ylabel('Displacement/m');

legend([fig1_comps.p1, fig1_comps.p2], 'windward', 'leeward');
legendX = .82; legendY = .87; legendWidth = 0.01; legendHeight = 0.01;
fig1_comps.legendPosition = [legendX, legendY, legendWidth, legendHeight];
% If you want the tightest box set width and height values very low matlab automatically sets the tightest box

%========================================================
% SET PLOT PROPERTIES
% Choices for COLORS can be found in ColorPalette.png
set(fig1_comps.p1, 'LineStyle', '-', 'LineWidth', 2, 'Color', PS.Blue4);
set(fig1_comps.p2, 'LineStyle', '-', 'LineWidth', 2, 'Color', PS.MyRed);

%========================================================
% INSTANTLY IMPROVE AESTHETICS-most important step
STANDARDIZE_FIGURE(fig1_comps);

%% 识别气动力参数
%% Identify aerodynamic parameters

% 结构参数（上游）
% Structural parameters （windward side）
D = 0.667; % deck depth
m = 80/L; % mass of the segment model per unit length
F0 = up_Fre_vibration; % Frequency without wind
omega0 = 2 * pi * F0; % Circular frequency without wind
rho = 1.225; % density of the air
% U = 6.21; % wind speed
Name = my_table.casename;
isexist = find(Name == fname);
if ~isempty(my_table.Windspeed(isexist))
    alertstr=(fname+"工况风速为："+num2str(my_table.Windspeed(isexist))+",是否需要修改？(0/不修改，1/修改)");
    decide=input(alertstr);
    if decide==0
        U=my_table.Windspeed(isexist);
    else
        alertstr="请输入该工况风速:";
        U= input(alertstr);
    end
else
    alertstr="请输入该工况风速:";
    U= input(alertstr);
end
% Zeta0 = 0.29/100; % damping ratio without wind

% 开始识别上游梁气动参数，并进行响应重构
[up_ex, up_frex] = ee(UP, 1 / fs); %经验包络法求瞬时频率和瞬时振幅 % empirical envelope method to find instantaneous frequency and instantaneous amplitude
up_ex = up_ex / D; %采用无量纲化 % Scale the envelope to dimensionless units
front_end = 0; end_end = 0;
slength = up_t(end);
up_frex = up_frex(fs * front_end + 1:fs * (slength - end_end)); %频率边界效应
up_ex = up_ex(fs * front_end + 1:fs * (slength - end_end)); %振幅边界效应
UP = UP(fs * front_end + 1:fs * (slength - end_end));
up_t = up_t(1:(slength - front_end - end_end) * fs);

% bex=polyfit(up_t,ex,7);%对振幅曲线进行拟合 已经改为 gompertz model
% ex=polyval(bex,up_t);
if fitdecide==1
[up_fitresult, up_gof] = createFit(up_t, up_ex);
figure
plot(up_t, up_ex)
hold on
up_ex = up_fitresult(up_t);
plot(up_t, up_ex)
xlabel("t / s")
ylabel("无量纲振幅")
legend("measured", "fitted")
title("无量纲振幅包络线拟合")
end

figure
plot(up_t, UP)
hold on
plot(up_t, up_ex * D, 'r')
xlabel("t / s")
ylabel("实际振幅")
legend("measured", "calculated")
title("原始振动数据与拟合包络线对比")

[up_psd_avg, up_f, up_psd_plot] = fft_transfer(fs, UP);

[up_psdmax, up_psdmaxseq] = max(up_psd_plot);
up_Fren_vibration_withwind = up_f(up_psdmaxseq); %计算系统振动频率
up_omegan_vibration_withwind = 2*pi*up_Fren_vibration_withwind; %计算系统振动频率对应的圆频率

% F0 = up_Fre_vibration; % Frequency without wind
omega0 = up_omega0_vibration; % Circular frequency without wind


up_omgx = up_frex * 2 * pi; %瞬时圆频率
up_bomgx = polyfit(up_ex, up_omgx, 4);
up_omgxeq = polyval(up_bomgx, up_ex); %瞬时圆频率多项式拟合
% figure; plot(up_ex, up_omgx, 'g'); hold on; plot(up_ex, up_omgxeq, 'r'); title('瞬时频率结果 计算结果(绿)+多项式拟合结果(红)+真实值(蓝)')

up_epsx = zeros(1, length(up_ex)); up_epsx = up_epsx'; %瞬时阻尼比

for i = 1:length(up_ex) - 1
    up_epsx(i) = log(up_ex(i) / up_ex(i + 1)) ./ up_omgx(i) * fs;
end

up_epsx(length(up_ex)) = up_epsx(length(up_ex) - 1);

decidevalue=input("是否需要限制上游0振幅阻尼比:(0/不需要，1/需要)");
if decidevalue==1
    up_epsx0=input("设置上游梁0振幅时阻尼比为：");
    up_ex=[0;up_ex];
    up_epsx=[up_epsx0;up_epsx];
end
up_bepsx = polyfit(up_ex, up_epsx, 4);
up_epsxeq = polyval(up_bepsx, up_ex); %瞬时阻尼比多项式拟合
figure; plot(up_ex * D, up_epsx, 'g'); hold on; plot(up_ex * D, up_epsxeq, 'r'); title('瞬时阻尼结果 计算结果(绿)+多项式拟合结果(红)')
% test=polyval(up_bepsx, up_amp)
up_amp = 0.0001:0.0001:1; up_amp = up_amp'; up_epsxeq = polyval(up_bepsx, up_amp);
figure
plot(up_amp * D, up_epsxeq)

for k1 = 1:length(up_bepsx)
    up_c(k1) = up_bepsx(end - k1 + 1); %拟合的多项式系数进行倒序排列
end

Zeta0 = Zeta0_up;
up_a1 = -2 * (up_c(1) - Zeta0) * omega0 * m / (rho * U * D);
up_a2 = -3 * up_c(2) * pi * omega0 * m / (2 * rho * U * D);
up_a3 = -8 * up_c(3) * omega0 * m / (rho * U * D);
up_a4 = -15 * up_c(4) * pi * omega0 * m / (4 * rho * U * D);
up_a5 = -16 * up_c(5) * omega0 * m / (rho * U * D);

up_a = [up_a1 up_a2 up_a3 up_a4 up_a5];
up_H4=(up_omega0_vibration^2-up_omegan_vibration_withwind^2)*m_up/rho/U^2;


Fre = up_Fre_vibration;
Mass = m;


dt = 1 / fs;
up_told = up_t;
up_tt = up_told;
fsup_a = [up_a1 up_a2 up_a3 up_a4 up_a5];
up_P = zeros(1, size(up_told, 1));
up_amp = 0.0001:0.0001:0.01; up_amp = up_amp';

up_epsall = -rho .* U .* D .* (up_a(1) + 4 .* up_a(2) .* up_amp ./ 3 ./ pi + up_a(3) .* up_amp.^2/4 + 8 .* up_a(4) .* up_amp.^3/15 / pi + up_a(5) .* up_amp.^4/8) ./ 2 ./ omega0 ./ m + Zeta0;
up_eps_aero = -rho .* U .* D .* (up_a(1) + 4 .* up_a(2) .* up_amp ./ 3 ./ pi + up_a(3) .* up_amp.^2/4 + 8 .* up_a(4) .* up_amp.^3/15 / pi + up_a(5) .* up_amp.^4/8) ./ 2 ./ omega0 ./ m;
figure
plot(up_amp * D, up_eps_aero)
figure
plot(up_amp * D, up_epsall)
hold on
plot(up_ex * D, up_epsx, 'g')

up_u0 = UP(1);
up_udot0 = (UP(2) - UP(1)) / dt;



% up_out = polynomial_NB_adstiff(Fre, Mass, Zeta0, rho, D, U, up_a, up_H4, up_told, up_P, up_u0, up_udot0);
% up_out = polynomial_NB(Fre, Mass, Zeta0, rho, D, U, up_a,  up_told, up_P, up_u0, up_udot0);

up_upperlimit=max(up_ex);%无量纲位移最大值
up_lowerlimit=min(up_ex);%无量纲位移最小值
up_out = polynomial_NB_adstiff_withlimit(Fre, Mass, Zeta0, rho, D, U, up_a, up_H4, up_told, up_P, up_u0, up_udot0,up_upperlimit,up_lowerlimit,up_Fre_vibration);

figure
plot(up_out(:, 1), up_out(:, 2))
hold on
plot(up_t, UP)
legend("calculated", "measured")
title("上游振动时程重构")


% [down_psd_avg_dltx, down_f_dltx, down_psd_plot_dltx] = fft_transfer(fs, UP);
% [down_psdmax_dltx, down_psdmaxseq_dltx] = max(down_psd_plot_dltx);
% down_Fre_vibration = down_f_dltx(down_psdmaxseq_dltx) %计算系统振动频率
% figure
% plot(down_f_dltx, down_psd_plot_dltx)
% hold on
% [down_psd_avg_dltx, down_f_dltx, down_psd_plot_dltx] = fft_transfer(fs, up_out(:, 2));
% [down_psdmax_dltx, down_psdmaxseq_dltx] = max(down_psd_plot_dltx);
% down_Fre_vibration = down_f_dltx(down_psdmaxseq_dltx) %计算系统振动频率
% plot(down_f_dltx, down_psd_plot_dltx)

clear Zeta0 Fre F0 omega0

% 结构参数（下游）
% Structural parameters （leeward side）
F0 = down_Fre_vibration; % Frequency without wind
omega0 = 2 * pi * F0; % Circular frequency without wind


% 开始识别下游梁气动参数，并进行响应重构
[down_ex, down_frex] = ee(DOWN, 1 / fs); %经验包络法求瞬时频率和瞬时振幅 % empirical envelope method to find instantaneous frequency and instantaneous amplitude
down_ex = down_ex / D; %采用无量纲化 % Scale the envelope to dimensionless units
front_end = 0; end_end = 0;
slength = down_t(end);
down_frex = down_frex(fs * front_end + 1:fs * (slength - end_end)); %频率边界效应
down_ex = down_ex(fs * front_end + 1:fs * (slength - end_end)); %振幅边界效应
DOWN = DOWN(fs * front_end + 1:fs * (slength - end_end));
down_t = down_t(1:(slength - front_end - end_end) * fs);


if fitdecide==1
    [down_fitresult, down_gof] = createFit(down_t, down_ex);
    figure
    plot(down_t, down_ex)
    hold on
    down_ex = down_fitresult(down_t);
    plot(down_t, down_ex)
    xlabel("t / s")
    ylabel("无量纲振幅")
    legend("measured", "fitted")
    title("无量纲振幅包络线拟合")
end

figure
plot(down_t, DOWN)
hold on
plot(down_t, down_ex * D, 'r')
xlabel("t / s")
ylabel("实际振幅")
legend("measured", "calculated")
title("原始振动数据与拟合包络线对比")


[down_psd_avg, down_f, down_psd_plot] = fft_transfer(fs, DOWN);

[down_psdmax, down_psdmaxseq] = max(down_psd_plot);

down_Fren_vibration_withwind = down_f(down_psdmaxseq); %计算系统振动频率
figure
plot(down_f, down_psd_plot)

down_omegan_vibration_withwind = 2*pi*down_Fren_vibration_withwind; %计算系统振动频率对应的圆频率
omega0 = down_omega0_vibration; % Circular frequency without wind

down_omgx = down_frex * 2 * pi; %瞬时圆频率
down_bomgx = polyfit(down_ex, down_omgx, 4);
down_omgxeq = polyval(down_bomgx, down_ex); %瞬时圆频率多项式拟合
% figure; plot(down_ex, down_omgx, 'g'); hold on; plot(down_ex, down_omgxeq, 'r'); title('瞬时频率结果 计算结果(绿)+多项式拟合结果(红)+真实值(蓝)')



down_epsx = zeros(1, length(down_ex)); down_epsx = down_epsx'; %瞬时阻尼比

for i = 1:length(down_ex) - 1
    down_epsx(i) = log(down_ex(i) / down_ex(i + 1)) ./ down_omgx(i) * fs;
end

down_epsx(length(down_ex)) = down_epsx(length(down_ex) - 1);

decidevalue=input("是否需要限制下游0振幅阻尼比:(0/不需要，1/需要)");
if decidevalue==1
    down_epsx0=input("设置下游梁0振幅时阻尼比为：");
    down_ex=[0;down_ex];
    down_epsx=[down_epsx0;down_epsx];
end

down_bepsx = polyfit(down_ex, down_epsx, 4);
down_epsxeq = polyval(down_bepsx, down_ex); %瞬时阻尼比多项式拟合
figure; plot(down_ex * D, down_epsx, 'g'); hold on; plot(down_ex * D, down_epsxeq, 'r'); title('瞬时阻尼结果 计算结果(绿)+多项式拟合结果(红)')

down_amp = 0.0001:0.0001:0.001; down_amp = down_amp'; down_epsxeq = polyval(down_bepsx, down_amp);
figure
plot(down_amp * D, down_epsxeq)

for k1 = 1:length(down_bepsx)
    down_c(k1) = down_bepsx(end - k1 + 1); %拟合的多项式系数进行倒序排列
end
Zeta0 = Zeta0_down;
down_a1 = -2 * (down_c(1) - Zeta0) * omega0 * m / (rho * U * D);
down_a2 = -3 * down_c(2) * pi * omega0 * m / (2 * rho * U * D);
down_a3 = -8 * down_c(3) * omega0 * m / (rho * U * D);
down_a4 = -15 * down_c(4) * pi * omega0 * m / (4 * rho * U * D);
down_a5 = -16 * down_c(5) * omega0 * m / (rho * U * D);

down_a = [down_a1 down_a2 down_a3 down_a4 down_a5];
down_H4=(down_omega0_vibration^2-down_omegan_vibration_withwind^2)*m_down/rho/U^2;

down_P = zeros(1, size(t, 1));
Fre = down_Fre_vibration;
Mass = m;
down_u0 = DOWN(1);
down_udot0 = (DOWN(2) - DOWN(1)) / (1 / fs);


dt = 1 / fs;
down_told = down_t;
down_tt = down_told;
fsdown_a = [down_a1 down_a2 down_a3 down_a4 down_a5];
down_P = zeros(1, size(down_told, 1));
down_amp =  0.0001:0.0001:0.01; down_amp = down_amp';


%NOTE:振幅之外插值的问题通过设定固定阻尼比解决。

down_epsall = -rho .* U .* D .* (down_a(1) + 4 .* down_a(2) .* down_amp ./ 3 ./ pi + down_a(3) .* down_amp.^2/4 + 8 .* down_a(4) .* down_amp.^3/15 / pi + down_a(5) .* down_amp.^4/8) ./ 2 ./ omega0 ./ m + Zeta0;

% down_bepsx_revised_fit = polyfit(down_amp, down_epsall_revised, 4);
% down_epsxeq_revised2 = polyval(down_bepsx_revised_fit, down_amp); %瞬时阻尼比多项式拟合
% figure; plot(down_amp * D, down_epsall_revised, 'g'); hold on; plot(down_amp * D, down_epsxeq_revised2, 'r'); title('瞬时阻尼结果 计算结果(绿)+多项式拟合结果(红)')
% down_amp = 0.0000:0.0001:0.05; down_amp = down_amp'; down_epsxeq = polyval(down_bepsx_revised_fit, down_amp);
figure
plot(down_amp * D, down_epsall)
hold on
plot(down_ex * D, down_epsx, 'g')


down_u0 = DOWN(1);
down_udot0 = (DOWN(2) - DOWN(1)) / dt;
% down_udot0 = 0;
down_tt = down_told;

% down_out = polynomial_NB_adstiff(Fre, Mass, Zeta0, rho, D, U, down_a,down_H4, down_told, down_P, down_u0, down_udot0);

down_upperlimit=max(down_ex);%无量纲位移最大值
down_lowerlimit=min(down_ex);%无量纲位移最小值

down_out = polynomial_NB_adstiff_withlimit(Fre, Mass, Zeta0, rho, D, U, down_a,down_H4, down_told, down_P, down_u0, down_udot0,down_upperlimit,down_lowerlimit,down_Fre_vibration);



figure
plot(down_out(:, 1), down_out(:, 2),'b')
hold on
plot(down_t, DOWN)
legend("calculated", "measured")
title("下游振动时程重构")
[down_psd_avg_dltx, down_f_dltx, down_psd_plot_dltx] = fft_transfer(fs, UP(end/2:end));
[down_psdmax_dltx, down_psdmaxseq_dltx] = max(down_psd_plot_dltx);
test1 = down_f_dltx(down_psdmaxseq_dltx) %计算系统振动频率
% figure
% plot(down_f_dltx, down_psd_plot_dltx)
% hold on
[down_psd_avg_dltx, down_f_dltx, down_psd_plot_dltx] = fft_transfer(fs, up_out(end/2:end, 2));
[down_psdmax_dltx, down_psdmaxseq_dltx] = max(down_psd_plot_dltx);
test2 = down_f_dltx(down_psdmaxseq_dltx) %计算系统振动频率
% plot(down_f_dltx, down_psd_plot_dltx)

[down_psd_avg_dltx, down_f_dltx, down_psd_plot_dltx] = fft_transfer(fs, DOWN(end/2:end));
[down_psdmax_dltx, down_psdmaxseq_dltx] = max(down_psd_plot_dltx);
test3 = down_f_dltx(down_psdmaxseq_dltx) %计算系统振动频率
% figure
% plot(down_f_dltx, down_psd_plot_dltx)
% hold on
[down_psd_avg_dltx, down_f_dltx, down_psd_plot_dltx] = fft_transfer(fs, down_out(end/2:end, 2));
[down_psdmax_dltx, down_psdmaxseq_dltx] = max(down_psd_plot_dltx);
test4 = down_f_dltx(down_psdmaxseq_dltx) %计算系统振动频率
% plot(down_f_dltx, down_psd_plot_dltx)


%%debug

% down_u0 = 1e-4;
% down_udot0 = 0;
% down_told=0:0.01:50;
% down_told=down_told';
% down_ex_upperlimit=100000000;
% down_ex_lowerlimit=0;
% down_P=zeros(1, size(down_told, 1));


% figure
% down_out2 = polynomial_NB(Fre, Mass, Zeta0, rho, D, U, down_a, down_told, down_P, down_u0, down_udot0);
% plot(down_out2(:, 1), down_out2(:, 2),'r')
% hold on
% down_out = polynomial_NB_withlimit(Fre, Mass, Zeta0, rho, D, U, down_a, down_ex_upperlimit, down_ex_lowerlimit, down_told, down_P, down_u0, down_udot0);
% plot(down_out(:, 1), down_out(:, 2),'b')
% 
% down_ex_lowerlimit=0.015;
% down_out = polynomial_NB_withlimit(Fre, Mass, Zeta0, rho, D, U, down_a, down_ex_upperlimit, down_ex_lowerlimit, down_told, down_P, down_u0, down_udot0);
% plot(down_out(:, 1), down_out(:, 2),'g')
% plot(down_t, DOWN)
% legend("calculated", "measured")
% title("下游振动时程重构")



% down_out = polynomial_NB(Fre, Mass, Zeta0, rho, D, U, down_a, down_told, down_P, down_u0, down_udot0);


% down_out2 = polynomial_NB(Fre, Mass, Zeta0, rho, D, U, down_a_new, down_told, down_P, down_u0, down_udot0);
% figure
% plot(down_out(:, 1), down_out(:, 2),'r')
% hold on
% plot(down_out2(:, 1), down_out2(:, 2),'b')
% legend("old", "new")

% close all
% 是否需要记录数据
if recorddata == 1
    T.casename = fname;
    T.case = casenumber;
    T.AOA = AOA;
    T.spacing = spacing;
    T.type = type;
    T.Rspeed = Rspeed;
    T.Windspeed = U;
    T.up_Fre_vibration=up_Fre_vibration;
    T.up_dltx_zeta0=Zeta0_up;
    T.up_Fren_vibration_withwind=up_Fren_vibration_withwind;
    T.up_ReducedFre=2*pi*up_Fre_vibration*D/U;
    T.down_Fre_vibration=down_Fre_vibration;
    T.down_dltx_zeta0=Zeta0_down;
    T.down_Fren_vibration_withwind=down_Fren_vibration_withwind;
    T.down_ReducedFre=2*pi*down_Fre_vibration*D/U;
    T.up_start = up_index_start;
    T.up_end = up_index_end;
    T.down_start = down_index_start;
    T.down_end = down_index_end;
    T.up_parameter_a1 = up_a(1);
    T.up_parameter_a2 = up_a(2);
    T.up_parameter_a3 = up_a(3);
    T.up_parameter_a4 = up_a(4);
    T.up_parameter_a5 = up_a(5);
    T.up_parameter_H4 = up_H4;
    T.down_parameter_a1 = down_a(1);
    T.down_parameter_a2 = down_a(2);
    T.down_parameter_a3 = down_a(3);
    T.down_parameter_a4 = down_a(4);
    T.down_parameter_a5 = down_a(5);
    T.down_parameter_H4 = down_H4;
    T.up_dltx_start=up_index_start_dltx;
    T.up_dltx_end=up_index_end_dltx;
    T.up_upperlimit=up_upperlimit;
    T.up_lowerlimit=up_lowerlimit;
    T.down_dltx_start=down_index_start_dltx;
    T.down_dltx_end=down_index_end_dltx;
    T.down_upperlimit=down_upperlimit;
    T.down_lowerlimit=down_lowerlimit;    

    Name = my_table.casename;
    isexist = find(Name == fname);

    if isempty(isexist)
        alertstr = "该工况尚未含有数据，是否添加该工况数据？（0/不添加；1/添加）";
        decidedata = input(alertstr);
        if decidedata == 1
            my_table_new=struct2table(T);
            my_table = [my_table; my_table_new];
            save SZTD110_logfile my_table;
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
            save SZTD110_logfile my_table;
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
