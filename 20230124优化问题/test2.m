
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

save salib.m result