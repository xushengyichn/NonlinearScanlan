%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: Shengyi xushengyichn@outlook.com
%Date: 2022-11-28 17:39:07
%LastEditors: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%LastEditTime: 2023-01-16 22:05:01
%FilePath: \NonlinearScanlan\20230116全阶模态对1TMD控制效果影响\data_analysis.m
%Description: 分析数据
%
%Copyright (c) 2022 by Shengyi xushengyichn@outlook.com, All Rights Reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 绘制分别考虑1-前5阶模态对TMD控制第一阶涡振效果的影响
clc; clear; close all;

% 读取数据
data = importdata('nmodes_onetmd_results_loc.mat');

loc=data(1:661,1);
dis=data(1:661,3);
figure
plot(loc,dis)

loc2=data(661*1+1:661*2,1);
dis2=data(661*1+1:661*2,3);
hold on
plot(loc2,dis2)

loc3=data(661*2+1:661*3,1);
dis3=data(661*2+1:661*3,3);
hold on
plot(loc3,dis3)

loc4=data(661*3+1:661*4,1);
dis4=data(661*3+1:661*4,3);
hold on
plot(loc4,dis4)

loc5=data(661*4+1:661*5,1);
dis5=data(661*4+1:661*5,3);
hold on
plot(loc5,dis5)

save nmodes_onetmd_results_loc_alldata.mat

%% 绘制考虑前8阶模态和考虑100阶的相对精确解的差异
clc
clear
close all

data = importdata('10modes_onetmd_results_loc.mat');
num=size(data,1)/9;
for k1 = 1:9
    str="dis"+num2str(data(1+num*(k1-1),2))+"=data(1+num*(k1-1):num*k1,3);";
    eval(str)
end
loc=data(1:num,1);

mode2effect=(dis2-dis1)./dis100;
mode3effect=(dis3-dis2)./dis100;
mode4effect=(dis4-dis3)./dis100;
mode5effect=(dis5-dis4)./dis100;
mode6effect=(dis6-dis5)./dis100;
mode7effect=(dis7-dis6)./dis100;
mode8effect=(dis8-dis7)./dis100;
mode100effect=(dis100-dis8)./dis100;

figure
plot(loc,mode2effect)
hold on
plot(loc,mode3effect)
% plot(loc,mode4effect)
% plot(loc,mode5effect)
% plot(loc,mode6effect)
% plot(loc,mode7effect)
% plot(loc,mode8effect)
% plot(loc,mode100effect)

figure
plot(loc,dis1)
hold on
plot(loc,dis2)
plot(loc,dis3)
% plot(loc,dis4)
% plot(loc,dis5)
% plot(loc,dis6)
% plot(loc,dis7)
% plot(loc,dis8)
% plot(loc,dis100)

modeinfo = load('modeinfo_all.mat');
nodegap=modeinfo.nodegap;
mode=modeinfo.mode_re;
for t1 = 1:length(loc)

    [~, index] = sort(abs(nodegap - loc(t1))); %查找与xTMD最接近的点的排序
    xResult = nodegap(index(1:2)); %获取最接近的两个点的x坐标
    mode2nodes = mode(index(1:2), 1:8); %获取两个点坐标的y值
    phi_result = interp1(xResult, mode2nodes, loc(t1), 'linear', 'extrap'); %插值以后任意点的振型
    modeTMD(t1, 1:8) = phi_result(1:8);
end
% figure
% plot(nodegap,mode(:,1))
% hold on
% plot(loc,modeTMD(:,1))

dis=[dis1 dis2 dis3 dis4 dis5 dis6 dis7 dis8 dis100];
modeeffect=[mode2effect mode3effect mode4effect mode5effect mode6effect mode7effect mode8effect];
Freq=modeinfo.Freq;

for k1 = 1:length(Freq)-1
    Freqeffect(k1,1)=(Freq(k1+1)-Freq(1))/Freq(1)*100;
end

% 对于60m处的点
pointseq=find(loc==60);
point_mode_effect=modeeffect(pointseq,:);
point_freq_effect=Freqeffect(1:7);
point_mode_shape=modeTMD(pointseq,:);

for k1 = 1:length(point_mode_shape)-1
   point_mode_shape_effect(k1,1)=abs((point_mode_shape(k1+1)-point_mode_shape(1))/point_mode_shape(1)*100);
end
figure
scatter3(point_mode_shape_effect,point_freq_effect,point_mode_effect)
xlabel("mode shape difference (%)")
ylabel("frequency differeence (%)")
zlabel("mode influence")


% 遍历所有点
figure
for k1 = 1:length(loc)
pointseq=find(loc==loc(k1));
point_mode_effect=modeeffect(pointseq,:);
point_freq_effect=Freqeffect(1:7);
point_mode_shape=modeTMD(pointseq,:);

for k2 = 1:length(point_mode_shape)-1
%    point_mode_shape_effect(k1,1)=((abs(point_mode_shape(k1+1))-abs(point_mode_shape(1)))/point_mode_shape(1)*100);
   point_mode_shape_effect(k2,1)=(abs(point_mode_shape(k2+1))/(max(abs(modeTMD(:,k2+1))))*100);
end
phi1=abs(point_mode_shape(1)/max(abs(modeTMD(:,1)))*ones(7,1))
scatter3(point_mode_shape_effect,point_freq_effect,point_mode_effect,50*ones(7,1),phi1)
hold on
end
xlabel("mode shape difference (%)")
ylabel("frequency differeence (%)")
zlabel("mode influence")

colorbar

% x=0:0.5:100;
% y=0:1:200;
% z=ones(201)*0.05;
% mesh(x,y,z, 'EdgeAlpha', 0.2)