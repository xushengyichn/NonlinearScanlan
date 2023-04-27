%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%Date: 2023-01-25 15:24:39
%LastEditors: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%LastEditTime: 2023-04-06 15:11:24
%FilePath: \TMD响应计算\Main.m
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

number_of_modes_to_control = [1]; % 需要控制的模态数 The number of modes to be controlled
number_of_modes_to_consider = 1; % 考虑的总模态数 The total number of modes considered
number_of_tmds = 1; % 考虑的总TMD数 The total number of TMDs
total_tmd_mass_ratio = 0.02; % 总质量比 The total mass ratio
mass_six_span = 10007779.7; % 深中通道非通航桥六跨连续梁质量 The mass of 6-span continuous beam of the non-navigational bridge of the Zhenzhong-Link
total_tmd_mass = total_tmd_mass_ratio * mass_six_span; % 总质量 The total mass
modal_damping_ratios = ones(1, number_of_modes_to_consider) * 0.003; % 模态阻尼 The damping of modes
t_length = 500; % 设定计算时间长度

% 输入初始参数
% Input initial parameters
TMDs_mass = [total_tmd_mass / number_of_tmds]; % 最后一个TMD的质量 The mass of the last TMD
TMDs_frequency = [0.83]; % TMD频率 The frequency of TMDs
TMDs_damping_ratio = [0.08]; % TMD阻尼 The damping of TMDs
TMDs_location = [55]; % TMD位置 The location of TMDs

[result] = a_0_main(number_of_modes_to_control, number_of_modes_to_consider, number_of_tmds, modal_damping_ratios, t_length, TMDs_mass, TMDs_frequency, TMDs_damping_ratio, TMDs_location);
t= result.t;
dis= result.dis;
figure
plot(t,dis)