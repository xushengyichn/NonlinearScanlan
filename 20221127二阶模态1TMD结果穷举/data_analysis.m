%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: Shengyi xushengyichn@outlook.com
%Date: 2022-11-28 17:39:07
%LastEditors: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%LastEditTime: 2023-01-16 11:27:11
%FilePath: \NonlinearScanlan\20221127二阶模态1TMD结果穷举\data_analysis.m
%Description: 分析数据
%
%Copyright (c) 2022 by Shengyi xushengyichn@outlook.com, All Rights Reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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

load modeinfo.mat
figure
plot(nodegap,mode(:,1))

% 针对276m位置号点研究其考虑1-5阶模态情况下响应的变化(位于nodegap中的第171行,同时位于dis中的277行)
caldata=zeros(5,5);
% caldata 为5行5列变量，分别为1-5阶模态下的振型、频率和响应，和每减少的一阶模态与考虑其模态对计算位移的影响,现在假设模型共5阶模态，因此第五阶模态假设为精确解
% 第五行表示与第一阶模态频率的差异
caldata(1,:)=mode(171,:);
caldata(2,:)=Freq;
caldata(3,:)=[dis(277),dis2(277),dis3(277),dis4(277),dis5(277)];
caldata(4,5)=(dis5(277)-dis4(277))/dis5(277);
caldata(4,4)=(dis4(277)-dis3(277))/dis5(277);
caldata(4,3)=(dis3(277)-dis2(277))/dis5(277);
caldata(4,2)=(dis2(277)-dis(277))/dis5(277);
caldata(4,1)=0;
caldata(5,1)=0;
caldata(5,2)=(Freq(2)-Freq(1))/Freq(1);
caldata(5,3)=(Freq(3)-Freq(2))/Freq(1);
caldata(5,4)=(Freq(4)-Freq(3))/Freq(1);
caldata(5,5)=(Freq(5)-Freq(4))/Freq(1);

figure
plotdata=caldata(:,2:end);
scatter3(abs(plotdata(1,:)),plotdata(5,:),plotdata(4,:))
xlabel("modal displacement")
ylabel("frequency difference")
zlabel("mode influence")