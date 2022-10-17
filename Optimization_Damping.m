%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%Date: 2022-10-14 11:41:51
%LastEditors: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%LastEditTime: 2022-10-17 20:26:57
%FilePath: \NonlinearScanlan\Optimization_Damping.m
%Description: 直接计算安装TMD后的阻尼比
%
%Copyright (c) 2022 by xushengyichn 54436848+xushengyichn@users.noreply.github.com, All Rights Reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc;clear;close all;
addpath("函数\")
calmodes=1;
mus=0.01/100:0.01/100:5/100;
Frequency=0.8:0.001:1.0;
Damping_ratio=0.05:0.001:0.25;
locationTMD=1001:1:1406;
[Ftmds_M,zetaTMDs_M,Massratio_M]=ndgrid(Frequency,Damping_ratio,locationTMD);
% 
% [Ftmds_M,zetaTMDs_M,Massratio_M]=ndgrid(Frequency,Damping_ratio,mus);
% variables = [Ftmds_M(:),zetaTMDs_M(:),Massratio_M(:)];
% 
% 
% numIterations = size(variables,1);
% res = zeros(numIterations, 1);


%% 单点计算
if 0
mass_six_span = 10007779.7;
mode_number=1;%气动力施加在第一阶模态
calmodes=1;
calmodes_all=[1];
numberofTMD=1;
mTMD=[mass_six_span*0.015];
zetaTMD=[0.1];
fTMD=[0.82];
locationTMD=[1179];

load opt_1tmd_5modes
mTMD=sol.mTMD;
zetaTMD=sol.zetaTMD;
fTMD=sol.fTMD;
locationTMD= sol.locationTMD;

% [result]=Compare_n_modes(mode_number,numberofTMD,mTMD,zetaTMD,fTMD,locationTMD,calmodes);
tic
[minDampingRatio,minDampingRatio_sys]=Optim_Damping_for_n_modes(mode_number,numberofTMD,mTMD,zetaTMD,fTMD,locationTMD,calmodes_all);
toc
disp(minDampingRatio)
disp(minDampingRatio_sys)
end

%% 穷举问题
mode_number=1;%气动力施加在第一阶模态
numberofTMD=1;
mass_six_span = 10007779.7;
zetaTMDs=0.05:0.001:0.25;
fTMDs=0.8:0.001:1.0;
mTMDs=0.015*mass_six_span;
locationTMDs=1001:1:1406;
[zetaTMDs_M,fTMDs_M,locationTMDs_M]=ndgrid(zetaTMDs,fTMDs,locationTMDs);
variables = [zetaTMDs_M(:),fTMDs_M(:),locationTMDs_M(:)];
numIterations=size(variables,1);
ppm = ParforProgressbar(size(variables,1));
parfor k1 = 1:numIterations
% for k1 = 1:numIterations
    zetaTMD=variables(k1,1);
    fTMD=variables(k1,2);
    locationTMD=variables(k1,3);
    [minDampingRatio]=Optim_Damping_for_n_modes(mode_number,numberofTMD,mTMDs,zetaTMD,fTMD,locationTMD,calmodes);
    res(k1)=minDampingRatio;
    pause(100/numIterations);
    ppm.increment();
end



%% 优化问题
mode_number=1;%气动力施加在第一阶模态
numberofTMD=1;
mass_six_span = 10007779.7;
prob = optimproblem('Description','Optimize the TMDs','ObjectiveSense','maximize'); %优化为最优阻尼比
zetaTMD =optimvar('zetaTMD',numberofTMD,'LowerBound',0.05,'UpperBound',0.25); %阻尼比
fTMD =optimvar('fTMD',numberofTMD,'LowerBound',0.8,'UpperBound',1.0); %频率
mTMD =optimvar('mTMD',numberofTMD,'LowerBound',0.001*mass_six_span,'UpperBound',0.015*mass_six_span); %质量
locationTMD = optimvar('locationTMD', numberofTMD,'Type', 'integer', 'LowerBound', 1001, 'UpperBound', 1406); %TMD所安装的节点

% options = optimoptions('simulannealbnd', 'Display', 'iter', 'PlotFcn', {'gaplotscorediversity', 'gaplotbestf', 'gaplotrange'});
options = optimoptions('ga', 'Display', 'iter', 'PlotFcn', {'gaplotscorediversity', 'gaplotbestf', 'gaplotrange'}, 'UseParallel', true);
[minDampingRatio,minDampingRatio_sys] = fcn2optimexpr(@Optim_Damping_for_n_modes,mode_number,numberofTMD,mTMD,zetaTMD,fTMD,locationTMD,[1]);
msum = sum(mTMD) == mass_six_span*0.015;

prob.Constraints.msum = msum;
% prob.Constraints.fTMD = fTMD ==0.82;
% prob.Constraints.locationTMD = locationTMD== 1179;

prob.Objective=minDampingRatio;
show(prob)

x0.zetaTMD = 0.1*ones(numberofTMD,1);
x0.fTMD = 0.82*ones(numberofTMD,1);
x0.mTMD = 0.015*mass_six_span*ones(numberofTMD,1);
x0.locationTMD = 1179*ones(numberofTMD,1);

[sol, optval] = solve(prob, x0, 'Options', options);

val = evaluate(prob.Objective,sol);
disp(val)
save opt_1tmd_5modes

