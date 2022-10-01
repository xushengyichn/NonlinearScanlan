%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%Date: 2022-09-14 10:30:54
%LastEditors: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%LastEditTime: 2022-09-14 12:19:35
%FilePath: \NonlinearScanlan\Compare_multicases.m
%Description: 试验数据与计算涡振振幅对比，多个工况，含TMD

%
%Copyright (c) 2022 by xushengyichn 54436848+xushengyichn@users.noreply.github.com, All Rights Reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;
% 试验数据路径
path(1) = "D:\资料存档\实验数据\2021年4月2日深中通道110m和60m连续梁涡振试验数据\测振\20210323\1.开槽间距6.7m\8.上下游吊臂上各装2个TMD，频率比≈1\+3\1.振幅曲线"; %台式机路径
path(2) = "D:\资料存档\实验数据\2021年4月2日深中通道110m和60m连续梁涡振试验数据\测振\20210323\1.开槽间距6.7m\2.有栏杆\+3\1.振幅曲线"; %台式机路径
path(3) = "D:\资料存档\实验数据\2021年4月2日深中通道110m和60m连续梁涡振试验数据\测振\20210323\1.开槽间距6.7m\2.有栏杆\+3\2.涡振期间气动力参数识别"; %台式机路径
path(4) = "C:\Users\shengyix\Documents\2021年4月2日深中通道110m和60m连续梁涡振试验数据\测振\20210323\1.开槽间距6.7m\2.有栏杆\+3\2.涡振期间气动力参数识别"; %NTNU笔记本路径
path(5) = "C:\Users\shengyix\Documents\2021年4月2日深中通道110m和60m连续梁涡振试验数据\测振\20210323\1.开槽间距6.7m\2.有栏杆\+3\dltx"; %NTNU笔记本路径 动力特性
path(6) = "C:\Users\shengyix\Documents\2021年4月2日深中通道110m和60m连续梁涡振试验数据\测振\20210323\1.开槽间距6.7m\8.上下游吊臂上各装2个TMD，频率比≈1\+3\1.振幅曲线";
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


ExpNames = [
    'SZTD-110-case8-22.3-2401';
    'SZTD-110-case8-22.3-2501';
    'SZTD-110-case8-22.3-2601';
    'SZTD-110-case8-22.3-2701';
    'SZTD-110-case8-22.3-2801';
    'SZTD-110-case8-22.3-2901';
    'SZTD-110-case8-22.3-3101';
    ]; %记录文件名
ExpNames2 = [
    'SZTD-110-case2-22.3-fasan-2401';
    'SZTD-110-case2-22.3-fasan-2501';
    'SZTD-110-case2-22.3-fasan-2601';
    'SZTD-110-case2-22.3-fasan-2701';
    'SZTD-110-case2-22.3-fasan-2801';
    'SZTD-110-case2-22.3-fasan-2901';
    'SZTD-110-case2-22.3-fasan-3101';
    ]; %记录文件名

girderindex=2;
if girderindex==1
    TMDsindex=[12 14 15 18];
    freshift=0.16;
elseif girderindex==2
    TMDsindex=[4 7 8 10];
    freshift=-0.128;
%     TMDsindex=[3 6 8 10];
%     freshift=-0.13;
end




for k1=1:size(ExpNames,1)    
% for k1=7
    ExpName=ExpNames(k1,:);
    ExpName2=ExpNames2(k1,:);

% 读取试验数据
[up_t, UP, down_t, DOWN] = ReadExpData(ExpName);

% PS = PLOT_STANDARDS();
% figure(1)
% fig1_comps.fig = gcf;
% hold on
% fig1_comps.p1 = plot(up_t, UP);
% fig1_comps.p2 = plot(down_t, DOWN);
% hold off

% 数值响应计算

[t_cal, dis_cal] = CalData_Polynomial_withTMD_singledegree(ExpName2, girderindex, TMDsindex, freshift);
% figure
% plot(t_cal,dis_cal)

% 计算与试验数据对比
if girderindex==1
    t_exp=up_t;
    dis_exp=UP;
elseif girderindex==2
    t_exp=down_t;
    dis_exp=DOWN;
end
figure
plot(t_exp,dis_exp)
hold on 
plot(t_cal,dis_cal)
% close all
% 计算与试验数据对比
expamp(k1)=std(dis_exp(round(end/2,0):end))*sqrt(2);
calamp(k1)=std(dis_cal(round(end/2,0):end))*sqrt(2);
end


result=table(ExpNames,expamp',calamp','VariableNames',{'ExpName','expamp','calamp'});
disp(result)
index=1:size(ExpNames,1);
figure
scatter(index,expamp)
hold on
scatter(index,calamp)
legend("Experiment","Calculation")


% disp("试验数据标准差："+expamp)
% disp("计算数据标准差："+calamp)
% close all