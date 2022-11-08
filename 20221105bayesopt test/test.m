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
fTMD1 = 0.8339;
% fTMD1 = optimizableVariable('fTMD1',[0.7 1.2],'Type','real');
xTMD1 = optimizableVariable('xTMD1',[0 660],'Type','real'); %TMD所安装的节点
% fTMD1 = optimizableVariable('fTMD1',[0 1],'Type','real');
% xTMD1 = optimizableVariable('xTMD1',[0 1],'Type','real'); %TMD所安装的节点


calmodes_all=1;
% mTMD=[mTMD1 mTMD2 mTMD3];
% zetaTMD=[zetaTMD1 zetaTMD2 zetaTMD3];
% fTMD=[fTMD1 fTMD2 fTMD3];
% xTMD=[xTMD1 xTMD2 xTMD3];

% fun=@(x)Optim_Damping_for_n_foces_n_modes_bayesopt2(mode_numbers,numberofTMD,x.fTMD1,x.xTMD1,calmodes_all,mu);
% 
% % fun2=@(X)xconstraint(X,mu,mass_six_span);
% results = bayesopt(fun,[fTMD1 xTMD1],'AcquisitionFunctionName','expected-improvement-per-second-plus','ExplorationRatio',0.5,'MaxObjectiveEvaluations',10000,'NumSeedPoints',1000);


fun=@(x)Optim_Damping_for_n_foces_n_modes_bayesopt(mode_numbers,numberofTMD,x.xTMD1,calmodes_all,mu);
results = bayesopt(fun,[xTMD1],'AcquisitionFunctionName','expected-improvement-per-second-plus','ExplorationRatio',0.5,'MaxObjectiveEvaluations',10000);
% results = bayesopt(fun,[xTMD1]);


