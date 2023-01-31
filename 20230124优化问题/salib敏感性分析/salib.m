%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%Date: 2023-01-27 17:39:55
%LastEditors: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%LastEditTime: 2023-01-30 10:09:57
%FilePath: \20230124优化问题\salib.m
%Description: 敏感性分析数据分析代码，用于后续读入python文件分析敏感性结果
%
%Copyright (c) 2023 by xushengyichn 54436848+xushengyichn@users.noreply.github.com, All Rights Reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

% 1. Load the data
param_values = load('param_values.txt');


% # # for i in range(10240):
% # #     mTMD = param_values[i,0]
% # #     fTMD = param_values[i,1]
% # #     TMD_damping_ratio = param_values[i,2]
% # #     pTMD = param_values[i,3]
% # #     test=eng.a_0_main(matlab.double([1,2,3,4,5,6]),10,1,matlab.double([0.003,0.003,0.003,0.003,0.003,0.003,0.003,0.003,0.003,0.003]),100,mTMD,fTMD,TMD_damping_ratio,pTMD)
% # #     Y[i]=test["dis_all_modes_sum"]
numIterations=size(param_values,1);
result=zeros(numIterations,1);
pauseTime = 60 / numIterations;
ppm = ParforProgressbar(numIterations,'showWorkerProgress',true,'progressBarUpdatePeriod',3,'title','my fancy title');

parfor k1 = 1:numIterations
    temp=param_values(k1,:);
    mTMD = temp(1);
    fTMD = temp(2);
    TMD_damping_ratio = temp(3);
    pTMD = temp(4);
    test=a_0_main([1 2 3 4 5 6], 10, 1, 0.003 * ones(10, 1), 150, mTMD, fTMD, TMD_damping_ratio, pTMD);
    result(k1)=test.dis_all_modes_sum;
    pause(pauseTime);
    % increment counter to track progress
    ppm.increment();

end

save salib.mat result
save salib.txt -ascii a