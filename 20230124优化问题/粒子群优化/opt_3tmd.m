%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%Date: 2023-01-30 16:08:01
%LastEditors: Shengyi xushengyichn@outlook.com
%LastEditTime: 2023-02-04 00:47:57
%FilePath: \NonlinearScanlan\20230124优化问题\粒子群优化\opt_3tmd.m
%Description: 优化3个TMD的参数，目标函数为最小振幅，目的是确认是否一个TMD可以控制两个模态
%
%Copyright (c) 2023 by ${git_name_email}, All Rights Reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%调用matlab自己的优化工具包计算

clc;clear all;close all;
addpath('..\')
addpath('..\..\函数\')



%% 优化问题(6模态5TMD)
if 1
total_tmd_mass_ratio = 0.02 % 总质量比 The total mass ratio
mass_six_span = 10007779.7 % 深中通道非通航桥六跨连续梁质量 The mass of 6-span continuous beam of the non-navigational bridge of the Zhenzhong-Link
total_tmd_mass = total_tmd_mass_ratio * mass_six_span % 总质量 The total mass
mTMD_single=total_tmd_mass/3 % 质量 The mass mTMD1
numberofTMD=3;
prob = optimproblem('Description','Optimize the TMD for mode 1','ObjectiveSense','maximize'); %优化为最优阻尼比
mTMD1 =optimvar('mTMD1',1,'LowerBound',mTMD_single*0.5,'UpperBound',mTMD_single*1.0); %频率
mTMD2 =optimvar('mTMD2',1,'LowerBound',mTMD_single*0.5,'UpperBound',mTMD_single*1.0); %频率
fTMD1 =optimvar('fTMD1',1,'LowerBound',0.7,'UpperBound',2); %频率
fTMD2 =optimvar('fTMD2',1,'LowerBound',0.7,'UpperBound',2); %频率
fTMD3 =optimvar('fTMD3',1,'LowerBound',0.7,'UpperBound',2); %频率
dTMD1 =optimvar('zetaTMD1',1,'LowerBound',0.05,'UpperBound',0.20); %阻尼比
dTMD2 =optimvar('zetaTMD2',1,'LowerBound',0.05,'UpperBound',0.20); %阻尼比
dTMD3 =optimvar('zetaTMD3',1,'LowerBound',0.05,'UpperBound',0.20); %阻尼比
xTMD1 = optimvar('xTMD1', 1,'LowerBound', 0, 'UpperBound', 660); %TMD所安装的节点
xTMD2 = optimvar('xTMD2', 1,'LowerBound', 0, 'UpperBound', 660); %TMD所安装的节点
xTMD3 = optimvar('xTMD3', 1,'LowerBound', 0, 'UpperBound', 660); %TMD所安装的节点

t_length=150;
number_of_modes_to_control=[1 2 3 4 5 6 ];
number_of_modes_to_consider=10;
number_of_tmds=3;
modal_damping_ratios=ones(1,number_of_modes_to_consider)*0.003;

% Define the output function



% options = optimoptions('particleswarm',Display='iter',PlotFcn='pswplotbestf');
% options = optimoptions('particleswarm',Display='iter',PlotFcn='pswplotbestf',UseParallel=true,OutputFcn=@storeTrial);
options = optimoptions('particleswarm',Display='iter',PlotFcn='pswplotbestf',OutputFcn=@storeTrial);
% options = optimoptions('ga', 'Display', 'iter', 'PlotFcn', {'gaplotscorediversity', 'gaplotbestf', 'gaplotrange'}, 'UseParallel', true);
result = fcn2optimexpr(@b_0_3_tmd,number_of_modes_to_control,number_of_modes_to_consider,number_of_tmds,modal_damping_ratios,t_length,mTMD1,mTMD2,fTMD1,fTMD2,fTMD3,dTMD1,dTMD2,dTMD3,xTMD1,xTMD2,xTMD3,total_tmd_mass);

% msum = sum(mTMD) == mass_six_span*0.015;
% prob.Constraints.msum = msum;

% prob.Constraints.fTMD = fTMD ==0.82;
% prob.Constraints.locationTMD = locationTMD== 1179;

prob.Objective=-result;
show(prob)



% x0.zetaTMD = 0.1*ones(numberofTMD,1);
% x0.fTMD = 0.82*ones(numberofTMD,1);
% x0.mTMD = 0.015*mass_six_span*ones(numberofTMD,1);
% x0.xTMD = 0.1*ones(numberofTMD,1);

[sol, optval] = solve(prob,'Solver','particleswarm','Options',options);

val = evaluate(prob.Objective,sol);
disp(val)

strmdoes=(num2str(number_of_modes_to_control));


str="save opt_pa_6mode_"+num2str(numberofTMD)+"TMD_"+strmdoes+"modes_xloc"+num2str(sol.xTMD)+".mat";
eval(str)
% save opt_1tmd_1modes
end


% Output function to store trial points
function stop = storeTrial(optimValues,state)
    x=optimValues.bestx;
    y=optimValues.bestfval;
    z=optimValues.iteration;
    x=[x y z];

    spmd
        filename = sprintf('opt_5tmd_%d.txt', labindex);
        if exist(filename, 'file')
            trial_points = dlmread(filename);
        else
            trial_points = x;
        end
        dlmwrite(filename, [trial_points; x], 'delimiter', ' ');
    end
    stop = false;
end
