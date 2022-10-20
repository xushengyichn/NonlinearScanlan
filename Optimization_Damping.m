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

% numberofTMD=1;
% mTMD=[mass_six_span*0.015];
% zetaTMD=[0.1];
% fTMD=[0.82];
% xTMD=56.2
load opt_1tmd_2modes
mTMD=sol.mTMD;
zetaTMD=sol.zetaTMD;
fTMD=sol.fTMD;
xTMD= sol.xTMD;


calmodes=2;
calmodes_all=[1];


[result]=Compare_n_modes(mode_number,numberofTMD,mTMD,zetaTMD,fTMD,xTMD,calmodes);
result.Mode
result.Mode_sys
% tic
% [minDampingRatio,minDampingRatio_sys]=Optim_Damping_for_n_modes(mode_number,numberofTMD,mTMD,zetaTMD,fTMD,xTMD,calmodes_all);
% toc
% disp(minDampingRatio)
% disp(minDampingRatio_sys)
end

%% 穷举问题(阻尼、频率、质量比)
if 0
mode_number=1;%气动力施加在第一阶模态
numberofTMD=1;
mass_six_span = 10007779.7;
zetaTMDs=0.05:0.001:0.25;
fTMDs=0.7:0.001:0.9;
muTMDs=0.003:0.003:0.03;
mTMDs=muTMDs*mass_six_span;
xTMDs=56.2;

[zetaTMDs_M,fTMDs_M,xTMDs_M]=ndgrid(zetaTMDs,fTMDs,mTMDs);
variables = [zetaTMDs_M(:),fTMDs_M(:),xTMDs_M(:)];
numIterations=size(variables,1);
ppm = ParforProgressbar(size(variables,1));
parfor k1 = 1:numIterations
% for k1 = 1:numIterations
    zetaTMD=variables(k1,1);
    fTMD=variables(k1,2);
    mTMD=variables(k1,3);
    xTMD=xTMDs;
    [minDampingRatio]=Optim_Damping_for_n_modes(mode_number,numberofTMD,mTMD,zetaTMD,fTMD,xTMD,calmodes);
    res(k1,1)=minDampingRatio;
%     disp(res)
    pause(100/numIterations);
    ppm.increment();
end

result=[variables res];

str="save traverse_"+num2str(numberofTMD)+"TMD_"+num2str(calmodes)+"modes_xloc"+num2str(xTMDs)+".mat result";
eval(str)

end

%% 穷举问题(阻尼、频率、位置)
if 0
mode_number=1;%气动力施加在第一阶模态
numberofTMD=1;
mass_six_span = 10007779.7;
zetaTMDs=0.05:0.001:0.25;
fTMDs=0.8:0.001:1.0;
mTMDs=0.015*mass_six_span;
xTMDs=0:1:606;
[zetaTMDs_M,fTMDs_M,xTMDs_M]=ndgrid(zetaTMDs,fTMDs,xTMDs);
variables = [zetaTMDs_M(:),fTMDs_M(:),xTMDs_M(:)];
numIterations=size(variables,1);
% ppm = ParforProgressbar(size(variables,1));
% parfor k1 = 1:numIterations
for k1 = 1:numIterations
    zetaTMD=variables(k1,1);
    fTMD=variables(k1,2);
    xTMD=variables(k1,3);
    [minDampingRatio]=Optim_Damping_for_n_modes(mode_number,numberofTMD,mTMDs,zetaTMD,fTMD,xTMD,calmodes);
    res(k1)=minDampingRatio;
%     disp(res)
%     pause(100/numIterations);
%     ppm.increment();
end
end



%% 优化问题(单模态一个TMD)
if 0
mode_number=1;%气动力施加在第一阶模态
numberofTMD=1;
mass_six_span = 10007779.7;
prob = optimproblem('Description','Optimize the TMDs','ObjectiveSense','maximize'); %优化为最优阻尼比
zetaTMD =optimvar('zetaTMD',numberofTMD,'LowerBound',0.05,'UpperBound',0.25); %阻尼比
fTMD =optimvar('fTMD',numberofTMD,'LowerBound',0.8,'UpperBound',1.0); %频率
mTMD =optimvar('mTMD',numberofTMD,'LowerBound',0.015*mass_six_span,'UpperBound',0.015*mass_six_span); %质量
xTMD = optimvar('xTMD', numberofTMD,'LowerBound', 56.2, 'UpperBound', 56.2); %TMD所安装的节点
calmodes_all=1:1:1;
% options = optimoptions('particleswarm',Display='iter',PlotFcn='pswplotbestf');
options = optimoptions('particleswarm',Display='iter',PlotFcn='pswplotbestf',UseParallel=true);
% options = optimoptions('ga', 'Display', 'iter', 'PlotFcn', {'gaplotscorediversity', 'gaplotbestf', 'gaplotrange'}, 'UseParallel', true);
[minDampingRatio,minDampingRatio_sys] = fcn2optimexpr(@Optim_Damping_for_n_modes,mode_number,numberofTMD,mTMD,zetaTMD,fTMD,xTMD,calmodes_all);

% msum = sum(mTMD) == mass_six_span*0.015;
% prob.Constraints.msum = msum;

% prob.Constraints.fTMD = fTMD ==0.82;
% prob.Constraints.locationTMD = locationTMD== 1179;

prob.Objective=minDampingRatio;
show(prob)

% x0.zetaTMD = 0.1*ones(numberofTMD,1);
% x0.fTMD = 0.82*ones(numberofTMD,1);
% x0.mTMD = 0.015*mass_six_span*ones(numberofTMD,1);
% x0.xTMD = 0.1*ones(numberofTMD,1);

[sol, optval] = solve(prob,'Solver','particleswarm','Options',options);

val = evaluate(prob.Objective,sol);
disp(val)

strmdoes=(num2str(calmodes_all));
strmdoes=strjoin(strsplit(strmdoes),'_');

str="save opt_"+num2str(numberofTMD)+"TMD_"+strmdoes+"modes_xloc"+num2str(sol.xTMD)+".mat";
eval(str)
% save opt_1tmd_1modes
end
%% 数据分析
if 0
str="load opt_"+num2str(numberofTMD)+"TMD_"+strmdoes+"modes_xloc"+num2str(sol.xTMD)+".mat";
eval(str)

mTMD=sol.mTMD;
zetaTMD=sol.zetaTMD;
fTMD=sol.fTMD;
xTMD= sol.xTMD;
[result]=Optim_Damping_for_n_modes(mode_number,numberofTMD,mTMD,zetaTMD,fTMD,xTMD,calmodes_all);
result
calmodes=2
[result]=Compare_n_modes(mode_number,numberofTMD,mTMD,zetaTMD,fTMD,xTMD,calmodes);
result.Mode
result.Mode_sys
end


%% 两个TMD
if 0
a=importdata("opt_1TMD_1_2modes_xloc56.2.mat")
b=importdata("opt_1TMD_1modes_xloc56.2.mat")
mode_number=1
numberofTMD=2
mTMD=[a.sol.mTMD b.sol.mTMD];
zetaTMD=[a.sol.zetaTMD b.sol.zetaTMD];
fTMD=[a.sol.fTMD b.sol.fTMD];
xTMD=[56.2 289.465]
calmodes=2
[result]=Compare_n_modes(mode_number,numberofTMD,mTMD,zetaTMD,fTMD,xTMD,calmodes);
result.Mode
result.Mode_sys


% test.fTMD=0.93
% test.mTMD=sol.mTMD;
% test.xTMD=56.2
% test.zetaTMD=0.22
% evaluate(prob.Objective,test);


% load damping_1tmd_5modes.mat
% 
% [~,seq]=max(result(:,4));
% opt_design=result(seq,:)
end

%% 单模态一个1个TMD（基于1 2 3 阶气动力）
calmodes_alls=[1 2 3];
mode_numbers=[1 2 3];

result=[];
for k1 = 1:3
mode_number=mode_numbers(k1);%气动力施加在第一阶模态
numberofTMD=1;
mass_six_span = 10007779.7;
prob = optimproblem('Description','Optimize the TMDs','ObjectiveSense','maximize'); %优化为最优阻尼比
zetaTMD =optimvar('zetaTMD',numberofTMD,'LowerBound',0.05,'UpperBound',0.25); %阻尼比
fTMD =optimvar('fTMD',numberofTMD,'LowerBound',0.8,'UpperBound',1.0); %频率
mTMD =optimvar('mTMD',numberofTMD,'LowerBound',0.015*mass_six_span,'UpperBound',0.015*mass_six_span); %质量
xTMD = optimvar('xTMD', numberofTMD,'LowerBound', 56.2, 'UpperBound', 56.2); %TMD所安装的节点
calmodes_all=calmodes_alls(k1);
options = optimoptions('particleswarm',Display='iter',PlotFcn='pswplotbestf',UseParallel=true);
[minDampingRatio,minDampingRatio_sys] = fcn2optimexpr(@Optim_Damping_for_n_modes,mode_number,numberofTMD,mTMD,zetaTMD,fTMD,xTMD,calmodes_all);


prob.Objective=minDampingRatio;
show(prob)


[sol, optval] = solve(prob,'Solver','particleswarm','Options',options);

val = evaluate(prob.Objective,sol);
disp(val)

strmdoes=(num2str(calmodes_all));
strmdoes=strjoin(strsplit(strmdoes),'_');

str="result.sol"+num2str(k1)+"=sol";
eval(str)
str="result.optval"+num2str(k1)+"=optval";
eval(str)
str="result.prob"+num2str(k1)+"=prob";
end