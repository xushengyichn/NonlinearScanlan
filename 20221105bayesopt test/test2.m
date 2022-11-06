%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%Date: 2022-11-05 23:44:42
%LastEditors: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%LastEditTime: 2022-11-05 23:53:42
%FilePath: \20221105bayesopt test\test2.m
%Description: 穷举所有工况，找到最优工况
%
%Copyright (c) 2022 by xushengyichn 54436848+xushengyichn@users.noreply.github.com, All Rights Reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear 
close all

addpath('../函数/')
mode_numbers=[1];%气动力施加的模态，输入n个数，表示分别计算n阶模态
numberofTMD=1;
mass_six_span = 10007779.7;
mu=0.02;
% mTMD1 = optimizableVariable('mTMD1',[0.01*mass_six_span 0.015*mass_six_span],'Type','real'); %质量
zetaTMD1 = 0.1;
% fTMD1 = 0.8339;
% fTMD1 = optimizableVariable('fTMD1',[0.7 1.2],'Type','real');;
% xTMD1 = optimizableVariable('xTMD1',[0 660],'Type','real'); %TMD所安装的节点
fTMD1=0.7:0.001:1.2;
xTMD1=0:1:660;

calmodes_all=1;

results=zeros(length(fTMD1),length(xTMD1));
[fTMD1_all,xTMD1_all]=ndgrid(fTMD1,xTMD1);
cases=[fTMD1_all(:),xTMD1_all(:)];

for k1 = 1:size(cases,1)
    fTMD1_temp=cases(k1,1);
    xTMD1_temp=cases(k1,2);
    [minDamping_allmodes,result]=Optim_Damping_for_n_foces_n_modes_bayesopt2(mode_numbers,numberofTMD,fTMD1_temp,xTMD1_temp,calmodes_all,mu);
    results(k1)=minDamping_allmodes;
end


% mTMD=[mTMD1 mTMD2 mTMD3];
% zetaTMD=[zetaTMD1 zetaTMD2 zetaTMD3];
% fTMD=[fTMD1 fTMD2 fTMD3];
% xTMD=[xTMD1 xTMD2 xTMD3];

% fun=@(x)Optim_Damping_for_n_foces_n_modes_bayesopt2(mode_numbers,numberofTMD,x.fTMD1,x.xTMD1,calmodes_all,mu);

% % fun2=@(X)xconstraint(X,mu,mass_six_span);
% results = bayesopt(fun,[fTMD1 xTMD1],'AcquisitionFunctionName','expected-improvement-per-second-plus','ExplorationRatio',0.5,'MaxObjectiveEvaluations',10000,'NumSeedPoints',1000);


% fun=@(x)Optim_Damping_for_n_foces_n_modes_bayesopt(mode_numbers,numberofTMD,x.xTMD1,calmodes_all,mu);
% results = bayesopt(fun,[xTMD1],'AcquisitionFunctionName','expected-improvement-per-second-plus','ExplorationRatio',0.5,'MaxObjectiveEvaluations',10000);
% % results = bayesopt(fun,[xTMD1]);


