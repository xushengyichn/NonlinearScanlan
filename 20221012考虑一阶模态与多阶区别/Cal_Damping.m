%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%Date: 2022-10-14 11:41:51
%LastEditors: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%LastEditTime: 2022-10-14 11:48:35
%FilePath: \NonlinearScanlan\20221012考虑一阶模态与多阶区别\Cal_Damping.m
%Description: 直接计算安装TMD后的阻尼比
%
%Copyright (c) 2022 by xushengyichn 54436848+xushengyichn@users.noreply.github.com, All Rights Reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc;clear;close all;
addpath("..\函数\")
calmodes=1;
mus=0.01/100:0.01/100:5/100;
Frequency=0.8:0.01:1.0;
Damping_ratio=0.05:0.01:0.15;

[Ftmds_M,zetaTMDs_M,Massratio_M]=ndgrid(Frequency,Damping_ratio,mus);
variables = [Ftmds_M(:),zetaTMDs_M(:),Massratio_M(:)];


numIterations = size(variables,1);
res = zeros(numIterations, 1);

ppm = ParforProgressbar(numIterations,'showWorkerProgress',true,'progressBarUpdatePeriod',3,'title','my fancy title'); 
pauseTime = 60/numIterations;
Sys_Fre=zeros(numIterations,2);
Sys_Damping=zeros(numIterations,2);
parfor k1 =1:size(variables,1)
% for k1 =100
    mu=variables(k1,3);
    zetaTMD=variables(k1,2);
    fTMD=variables(k1,1);
    locationTMD=1036;    
    [Mode]=Damping_n_modes(mu,zetaTMD,fTMD,locationTMD,calmodes);
    Sys_Fre(k1,:)=Mode.Frequency';
    Sys_Damping(k1,:)=Mode.("Damping ratio")';
    pause(pauseTime);
    % increment counter to track progress
    ppm.increment();
end

result=[variables Sys_Fre Sys_Damping];