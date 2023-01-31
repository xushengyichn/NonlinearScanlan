%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%Date: 2023-01-25 15:24:39
%LastEditors: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%LastEditTime: 2023-01-31 21:58:36
%FilePath: \20230124优化问题\traverse_mode_1to6.m
%Description: 本代码通过穷举法遍历所有可能的参数组合，得到最优的参数组合
%
%Copyright (c) 2023 by xushengyichn 54436848+xushengyichn@users.noreply.github.com, All Rights Reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear
close all
addpath("../函数/")

for t1 = 1:6
    % 输入预先设定的参数
    % Preset parameters
    total_tmd_mass_ratio = 0.02; % 总质量比 The total mass ratio
    mass_six_span = 10007779.7; % 深中通道非通航桥六跨连续梁质量 The mass of 6-span continuous beam of the non-navigational bridge of the Zhenzhong-Link
    total_tmd_mass = total_tmd_mass_ratio * mass_six_span; % 总质量 The total mass

    % 输入初始参数
    % Input initial parameters
    t_length = 150; % 设定计算时间长度
    number_of_tmds = 1; % 考虑的总TMD数 The total number of TMDs
    number_of_modes_to_consider = 10; % 考虑的总模态数 The total number of modes considered
    modal_damping_ratios = ones(1, number_of_modes_to_consider) * 0.003; % 模态阻尼 The damping of modes
    number_of_modes_to_control = [t1]; % 需要控制的模态数 The number of modes to be controlled
    TMDs_mass = [total_tmd_mass / number_of_tmds / 6]; % TMD质量 The mass of TMDs
    % % TMDs_mass = [TMDs_mass, total_tmd_mass - sum(TMDs_mass)]; % 最后一个TMD的质量 The mass of the last TMD
    % TMDs_frequency = [0.8, 1, 1.2]; % TMD频率 The frequency of TMDs
    % TMDs_damping_ratio = [0.08, 0.09, 0.07]; % TMD阻尼 The damping of TMDs
    % TMDs_location = [55, 155, 255]; % TMD位置 The location of TMDs

    f_lim = [0.8 1.0; 0.85 1.05; 0.95 1.15; 1.15 1.4; 1.45 1.6; 1.6 1.8];

    mtmd1 = total_tmd_mass / number_of_tmds / 6;
    ftmd1 = f_lim(t1, 1):0.01:f_lim(t1, 2);
    % ftmd1 = 0.8:0.01:1.0;
    dtmd1 = 0.05:0.05:0.20;
    xtmd1 = 0:3:660;

    [Ftmd1 Dtmd1 Xtmd1] = ndgrid(ftmd1, dtmd1, xtmd1);
    variables = [Ftmd1(:) Dtmd1(:) Xtmd1(:)];
    size(variables, 1) / 3600
    numIterations = size(variables, 1);

    collectdata = zeros(numIterations, 1);
    pauseTime = 60 / numIterations;
    ppm = ParforProgressbar(numIterations, 'showWorkerProgress', true, 'progressBarUpdatePeriod', 3, 'title', 'my fancy title');

    parfor k1 = 1:numIterations
        TMDs_mass = [mtmd1];
        TMDs_frequency = [variables(k1, 1)];
        TMDs_damping_ratio = [variables(k1, 2)];
        TMDs_location = [variables(k1, 3)];
        [result] = a_0_main(number_of_modes_to_control, number_of_modes_to_consider, number_of_tmds, modal_damping_ratios, t_length, TMDs_mass, TMDs_frequency, TMDs_damping_ratio, TMDs_location);
        collectdata(k1) = result.dis_all_modes_sum;
        pause(pauseTime);
        % increment counter to track progress
        ppm.increment();
    end

    del(ppm)
    str=['traverse_mode_',num2str(t1),'.mat'];
    save(str,'collectdata','variables')
end
