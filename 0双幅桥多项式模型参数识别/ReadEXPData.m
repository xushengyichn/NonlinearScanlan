%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: Shengyi Xu xushengyichn@outlook.com
%Date: 2022-07-01 16:26:47
%LastEditors: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%LastEditTime: 2022-08-29 19:43:10
%FilePath: \twindeck_ID\ReadEXPData.m
%Description: 程序功能是读取任意试验时程并画出图像The function of the program is to read any test data and draw an image
%
%Copyright (c) 2022 by Shengyi Xu xushengyichn@outlook.com, All Rights Reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;

%% 各种参数设置
%NTNU笔记本路径
addpath(genpath("C:\Users\shengyix\OneDrive\NAS云同步\Drive\0研究生\有用的代码\HHT-Tutorial-master")) %经验包络法求瞬时频率和振幅必备
addpath("C:\Users\shengyix\OneDrive\NAS云同步\Drive\0研究生\有用的代码\Matlab_PlottingTemplates") %绘图工具
%台式机路径
addpath(genpath("C:\Users\xushe\OneDrive\NAS云同步\Drive\0研究生\有用的代码\HHT-Tutorial-master"))
addpath("C:\Users\xushe\OneDrive\NAS云同步\Drive\0研究生\有用的代码\Matlab_PlottingTemplates")

%% 添加试验数据路径
path(1) = "C:\Users\shengyix\Documents\2021年4月2日深中通道110m和60m连续梁涡振试验数据\测振\20210323\1.开槽间距6.7m\2.有栏杆\+3\2.涡振期间气动力参数识别"; %NTNU笔记本路径
% path(2) = "D:\资料存档\实验数据\2021年4月2日深中通道110m和60m连续梁涡振试验数据\测振\20210323\1.开槽间距6.7m\2.有栏杆\+3\dltx"; %台式机路径
path(2) = "D:\资料存档\实验数据\2021年4月2日深中通道110m和60m连续梁涡振试验数据\测振\20210323\1.开槽间距6.7m\8.上下游吊臂上各装2个TMD，频率比≈1\+3\1.振幅曲线"; %台式机路径
% path(2) = "D:\资料存档\实验数据\2021年4月2日深中通道110m和60m连续梁涡振试验数据\测振\20210323\TMD标定\1"; %台式机路径
% path(3) = "D:\资料存档\实验数据\2021年4月2日深中通道110m和60m连续梁涡振试验数据\测振\20210323\1.开槽间距6.7m\8.上下游吊臂上各装2个TMD，频率比≈1\+3\1.振幅曲线";

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


%% 提取工况设置
% fname = ['SZTD-110-case2-22.3-fasan-2401']; %记录文件名
% fname = ['SZTD-110-case8-22.3-2401']; %记录文件名
% fname = ['SZTD-110-TMD1-8']; %记录文件名
% fname = ['SZTD-110-case2-22.3-dltx4']; %记录文件名
fname = ['SZTD-110-case8-22.3-3101']; %记录文件名

chanum = 8; %记录通道数
% filename = strsplit(fname, '-');
% casenumber = cell2mat(filename(3));
% spacing = str2double(cell2mat(filename(4)));
% type = cell2mat(filename(5));
% Rspeed = cell2mat(filename(6));

% casenumber = str2double(casenumber(5:end));
% Rspeed = Rspeed(1:end - 1);
% Rspeed = str2double(Rspeed);

% if strcmp(type, 'fasan')
%     type = 1;
% else

%     if strcmp(type, 'shuaijian')
%         type = 2;
%     else
%         error("实验数据文件名可能有误，请确认信号状态为衰减或发散，程序终止")
%     end

% end

% % 判断工况风攻角
% if casenumber <= 17 || and(casenumber >= 32, casenumber <= 35) || and(casenumber >= 40, casenumber <= 16) || casenumber == 55
%     AOA = 3;
% else

%     if and(casenumber >= 18, casenumber <= 25) || and(casenumber >= 36, casenumber <= 37) || and(casenumber >= 47, casenumber <= 50) || casenumber == 54
%         AOA = 0;
%     else

%         if and(casenumber >= 26, casenumber <= 31) || and(casenumber >= 38, casenumber <= 39) || and(casenumber >= 51, casenumber <= 53) || and(casenumber >= 56, casenumber <= 57)
%             AOA = -3;
%         end

%     end

% end

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

% save test.mat up_t UP

close all
% 开始计算下游振动时程重构
disp("基于TMD将振幅放大四倍，由于仅使用了一个激光位移计")
DOWN=DOWN*4;
[down_ex, down_frex] = ee(DOWN, 1 / fs); %经验包络法求瞬时频率和瞬时振幅 % empirical envelope method to find instantaneous frequency and instantaneous amplitude
front_end = 4; end_end = 2;
slength = down_t(end);
down_frex = down_frex(fs * front_end + 1:fs * (slength - end_end)); %频率边界效应
down_ex = down_ex(fs * front_end + 1:fs * (slength - end_end)); %振幅边界效应
DOWN = DOWN(fs * front_end + 1:fs * (slength - end_end));
down_t = down_t(1:(slength - front_end - end_end) * fs);
DOWN=DOWN-mean(DOWN);


figure
plot(down_t, DOWN)
[down_ex, down_frex] = ee(DOWN, 1 / fs); %经验包络法求瞬时频率和瞬时振幅 % empirical envelope method to find instantaneous frequency and instantaneous amplitude
hold on 


fitamp=polyfit(down_t,down_ex,4);
down_ex=polyval(fitamp,down_t);%Ë²Ê±Ô²ÆµÂÊ¶àÏîÊ½ÄâºÏ

plot(down_t, down_ex)

[down_psd_avg, down_f, down_psd_plot] = fft_transfer(fs, DOWN);
figure 
plot(down_f, down_psd_plot)
[down_psdmax, down_psdmaxseq] = max(down_psd_plot);
down_Fre_vibration = down_f(down_psdmaxseq); %计算系统振动频率

down_omgx = down_frex * 2 * pi; %瞬时圆频率
down_bomgx = polyfit(down_ex, down_omgx, 4);
down_omgxeq = polyval(down_bomgx, down_ex); %瞬时圆频率多项式拟合
% figure; plot(down_ex, down_omgx, 'g'); hold on; plot(down_ex, down_omgxeq, 'r'); title('瞬时频率结果 计算结果(绿)+多项式拟合结果(红)+真实值(蓝)')
figure; plot(down_ex, down_frex, 'g'); hold on; plot(down_ex, down_omgxeq/2/pi, 'r'); title('瞬时频率结果 计算结果(绿)+多项式拟合结果(红)+真实值(蓝)')

down_epsx = zeros(1, length(down_ex)); down_epsx = down_epsx'; %瞬时阻尼比

for i = 1:length(down_ex) - 1
    down_epsx(i) = log(down_ex(i) / down_ex(i + 1)) ./ down_omgx(i) * fs;
end

down_epsx(length(down_ex)) = down_epsx(length(down_ex) - 1);
down_bepsx = polyfit(down_ex, down_epsx, 4);
down_epsxeq = polyval(down_bepsx, down_ex); %瞬时阻尼比多项式拟合
figure; plot(down_ex, down_epsx, 'g'); hold on; plot(down_ex, down_epsxeq, 'r'); title('瞬时阻尼结果 计算结果(绿)+多项式拟合结果(红)')

down_t=down_t(end*0.08:end*0.25);
down_t=down_t-down_t(1);
DOWN=DOWN(end*0.5:end);

[down_psd_avg, down_f, down_psd_plot] = fft_transfer(fs, DOWN);
figure 
plot(down_f, down_psd_plot)
[down_psdmax, down_psdmaxseq] = max(down_psd_plot);
down_Fre_vibration = down_f(down_psdmaxseq); %计算系统振动频率
xlim([4, 6])
% t=0:1/fs:4;
% 
% dis=-0.05*cos(5.27*2*pi*(t-0.05));
% figure
% plot(down_t,DOWN)
% hold on
% plot(t,dis)
clear DOWN

%========================================================

Files=dir(strcat(path(2),'\*.sts'));  %列出指定目录下后缀为.sts的文件;
N=length(Files);
for k1 = 1:N
    DAT(k1,:)=Files(k1).name;
end
clear k1 N
DAT=string(DAT);
for i=1:length(DAT)
    A(i,:) = strsplit(DAT(i),'#');
end
B=A(:,1);
C=unique(B);
for k1 = 1:length(C)
    Csplit(k1,:)=strsplit(C(k1,1),'-');
end
Rspeed=Csplit(:,5)
Rspeed=(double(Rspeed)-1)/10;

for k1 = 1:length(C)
    fname = C(k1);
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

    DOWN=DOWN*4;
    [down_psd_avg, down_f, down_psd_plot] = fft_transfer(fs, DOWN);
    [a,b]=max(down_psd_plot);
    Fre(k1)=down_f(b);
    figure
    plot(down_t,DOWN)
end

figure
plot(Rspeed*0.027,Fre)
for k1 = 1:length(Fre)
    text(Rspeed(k1)*0.027,Fre(k1),num2str(round(Fre(k1),1)))
end
xlabel("wind speed")
ylabel("frequency")