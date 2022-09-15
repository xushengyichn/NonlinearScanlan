%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%Date: 2022-08-29 12:37:58
%LastEditors: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%LastEditTime: 2022-09-12 21:02:17
%FilePath: \NonlinearScanlan\Compare.m
%Description: 计算振幅与试验振幅对比，仅为单个工况对比，后续调整为批量对比
%
%Copyright (c) 2022 by xushengyichn 54436848+xushengyichn@users.noreply.github.com, All Rights Reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc; clear; close all;


test_Polynomial_withTMD_singledegree
clearvars -except  out Result2

%% 添加试验数据路径
path(1) = "D:\资料存档\实验数据\2021年4月2日深中通道110m和60m连续梁涡振试验数据\测振\20210323\1.开槽间距6.7m\8.上下游吊臂上各装2个TMD，频率比≈1\+3\1.振幅曲线"; %台式机路径
path(2) = "D:\资料存档\实验数据\2021年4月2日深中通道110m和60m连续梁涡振试验数据\测振\20210323\1.开槽间距6.7m\2.有栏杆\+3\1.振幅曲线"; %台式机路径
path(3) = "D:\资料存档\实验数据\2021年4月2日深中通道110m和60m连续梁涡振试验数据\测振\20210323\1.开槽间距6.7m\2.有栏杆\+3\2.涡振期间气动力参数识别"; %台式机路径


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
fname = ['SZTD-110-case8-22.3-2401']; %记录文件名
% fname = ['SZTD-110-case2-22.3-fasan-3101']; %记录文件名
chanum = 8; %记录通道数
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
% fig1_comps.p2 = plot(down_t, DOWN);
hold off

close all
figure
plot(up_t, UP);
hold on
plot(out(:, 1), out(:, 2))
legend('experiment','calculate')
disp(Result2)



caldata=out(round(end/2,0):end,2);
[psd_avg, f2, psd_plot2] = fft_transfer(256,caldata);
disp("计算主结构响应最大值"+num2str(std(out(end/2:end, 2)*sqrt(2))))
[a,b]=max(psd_plot2);
disp("计算振动频率为："+num2str(f2(b)))
% clearvars a b 
[psd_avg, f3, psd_plot3] = fft_transfer(1/(up_t(2)-up_t(1)),UP);
disp("试验测得响应均最大值"+num2str(std(UP(end/2:end/3*2))*sqrt(2)))
[c,d]=max(psd_plot3);
disp("试验测得振动频率为："+num2str(f3(d)))
disp(std(UP(end/2:end/3*2))*sqrt(2))
disp(std(out(end/2:end, 2)*sqrt(2)))
disp(f3(d))
disp(f2(b))
data(1)=std(UP(end/2:end/3*2))*sqrt(2);
data(2)=std(out(end/2:end, 2)*sqrt(2));
data(3)=f3(d);
data(4)=f2(b);