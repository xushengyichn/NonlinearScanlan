%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%Date: 2023-01-25 15:24:39
%LastEditors: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%LastEditTime: 2023-01-25 16:46:01
%FilePath: \20230124优化问题\Dataset_simple_question.m
%Description: 本脚本功能为输入一组TMD参数以及气动力控制参数，获得多阶模态涡振的响应
%
%Copyright (c) 2023 by xushengyichn 54436848+xushengyichn@users.noreply.github.com, All Rights Reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear
close all
addpath("../函数/")
% 输入预先设定的参数
% Preset parameters
total_tmd_mass_ratio = 0.02; % 总质量比 The total mass ratio
mass_six_span = 10007779.7; % 深中通道非通航桥六跨连续梁质量 The mass of 6-span continuous beam of the non-navigational bridge of the Zhenzhong-Link
total_tmd_mass = total_tmd_mass_ratio * mass_six_span; % 总质量 The total mass

% 输入初始参数
% Input initial parameters
t_length = 100; % 设定计算时间长度
number_of_tmds = 2; % 考虑的总TMD数 The total number of TMDs
number_of_modes_to_consider = 10; % 考虑的总模态数 The total number of modes considered
modal_damping_ratios = ones(1, number_of_modes_to_consider) * 0.003; % 模态阻尼 The damping of modes
number_of_modes_to_control = [1, 2, 3, 4, 5, 6]; % 需要控制的模态数 The number of modes to be controlled
TMDs_mass = [total_tmd_mass / number_of_tmds, total_tmd_mass / number_of_tmds]; % TMD质量 The mass of TMDs
TMDs_mass = [TMDs_mass, total_tmd_mass - sum(TMDs_mass)]; % 最后一个TMD的质量 The mass of the last TMD
TMDs_frequency = [0.8, 1, 1.2]; % TMD频率 The frequency of TMDs
TMDs_damping_ratio = [0.08, 0.09, 0.07]; % TMD阻尼 The damping of TMDs
TMDs_location = [55, 155, 255]; % TMD位置 The location of TMDs

mtmd1=1e4:10000:total_tmd_mass/2;
ftmd1=0.80:0.02:0.95;
dtmd1=0.08:0.02:0.1;
xtmd1=0:20:330;
ftmd2=0.80:0.02:0.95;
dtmd2=0.08:0.02:0.1;
xtmd2=0:20:330;

[Mtmd1 Ftmd1 Dtmd1 Xtmd1 Ftmd2 Dtmd2 Xtmd2]=ndgrid(mtmd1,ftmd1,dtmd1,xtmd1,ftmd2,dtmd2,xtmd2);
variables=[Mtmd1(:) Ftmd1(:) Dtmd1(:) Xtmd1(:) Ftmd2(:) Dtmd2(:) Xtmd2(:)];
size(variables,1)/3600/24*10/20
numIterations = size(variables, 1);

collectdata=zeros(numIterations,1);

ppm = ParforProgressbar(numIterations, 'showWorkerProgress', true, 'progressBarUpdatePeriod', 3, 'title', 'my fancy title');

for k1 = 1:numIterations
    TMDs_mass=[variables(k1,1),total_tmd_mass - sum(TMDs_mass)];
    TMDs_frequency=[variables(k1,2),variables(k1,5)];
    TMDs_damping_ratio=[variables(k1,3),variables(k1,6)];
    TMDs_location=[variables(k1,4),variables(k1,7)];
    [result] = a_0_main(number_of_modes_to_control, number_of_modes_to_consider, number_of_tmds, modal_damping_ratios, t_length, TMDs_mass, TMDs_frequency, TMDs_damping_ratio, TMDs_location);
    collectdata(k1)=dis_all_modes_sum
    pause(pauseTime);
        % increment counter to track progress
    ppm.increment();
end
