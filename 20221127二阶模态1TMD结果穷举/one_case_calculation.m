%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%Date: 2023-03-11 16:47:50
%LastEditors: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%LastEditTime: 2023-03-11 16:53:25
%FilePath: \NonlinearScanlan\20221127二阶模态1TMD结果穷举\one_case_calculation.m
%Description: 验证二阶模态下，安装一个TMD是否会对三阶模态都提高阻尼比
%
%Copyright (c) 2023 by ${git_name_email}, All Rights Reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc; clear; close all % 清除记录
addpath("../函数/")

numberofTMD = 1; % 所需要计算的TMD的数量.

savedata = 0;

nTMD = numberofTMD;
% nTMD = 2; %TMD数量
modeinfo = load('modeinfo.mat');
% 设计第一个TMD的参数
% TMD1按照规范设计
fs = modeinfo.Freq(1);
mode1 = modeinfo.eig_vec(:, 1);
mode2 = modeinfo.eig_vec(:, 2);
mu = 0.02;
mTMD1 = mu / (max(mode1)^2);
fTMD1 = 1 / (1 + mu) * fs;
zetaTMD1 = sqrt(3 * mu / 8 / (1 + mu));

xTMD2_all = 275;
fTMD2_all = 0.833;
zetaTMD2_all = 0.043;

[XTMD2_all, FTMD2_all, ZetaTMD2_all] = ndgrid(xTMD2_all, fTMD2_all, zetaTMD2_all);
variables = [XTMD2_all(:), FTMD2_all(:), ZetaTMD2_all(:)];

twomode_onetmd_dis = zeros(size(variables, 1), 4); %第四列为是否收敛标志
numIterations = size(variables, 1);

% ppm = ParforProgressbar(numIterations, 'showWorkerProgress', true, 'progressBarUpdatePeriod', 3, 'title', 'my fancy title');

pauseTime = 60 / numIterations;


% parfor k1 = 1:size(variables, 1)
for k1 = 1:1

    mass_six_span = 10007779.7;
    mTMD = [mTMD1]; %该代码为基准tmd的
    Ftmd = [variables(k1, 2)];
    zetaTMD = [variables(k1, 3)];
    xTMD = [variables(k1, 1)];

    omegaTMD = 2 * pi * Ftmd;

    mode_numbers = 2;
    ifcalmode = 3;
    h = 0.01;
    t_length = 150;

    %% 计算不安装TMD情况下各阶模态各点最大位移
    nodeondeck = importdata('nodeondeck.txt');
    KMmapping = importmappingmatrix('KMatrix.mapping');
    t = 0:h:t_length; % Time

    mode_number = mode_numbers;
    iter = 1;
    flag_convergence = 0;
        calmodes = mode_number; %考虑模态数 Consider the number of modes
        nModes = calmodes;

        MM_eq = modeinfo.MM_eq(1:calmodes, 1:calmodes);
        KK_eq = modeinfo.KK_eq(1:calmodes, 1:calmodes);
        eig_val = modeinfo.eig_val(1:calmodes, 1:calmodes);
        eig_vec = modeinfo.eig_vec(:, 1:calmodes);

        % 计算桥面节点的振型向量

        UYNode = sort(KMmapping.Node(KMmapping.DOF == 'UY'));
        mode = zeros(length(nodeondeck), nModes);
        mode_vec = eig_vec;

        for k3 = 1:length(nodeondeck)
            position_index = KMmapping.MatrixEqn(find(and(KMmapping.Node == nodeondeck(k3), KMmapping.DOF == 'UY')));

            if isempty(position_index)
                mode(k3, :) = zeros(1, nModes);
            else
                mode(k3, :) = mode_vec(position_index, :);
            end

        end

        mode_re = zeros(size(mode, 1), size(mode, 2));

        for k2 = 1:nModes
            mode_re(:, k2) = mode(:, k2) ./ max(abs(mode(:, k2)));
        end

        [phideckmax, location] = max(mode); %计算桥面每个模态的最大振型值
        nodegap = importdata('nodegap.txt');
        modes = [nodegap nodeondeck mode mode_re];
        %      if length(mTMD)==1
        %     for k1 = 1:length(mTMD)
        %         modeTMD(k1,1:mode_number)=mode(nodeTMD(k1)==modes(:,2),:);
        %     end

        modeTMD = zeros(length(mTMD), mode_number);

        for t1 = 1:length(mTMD)

            [~, index] = sort(abs(nodegap - xTMD)); %查找与xTMD最接近的点的排序
            xResult = nodegap(index(1:2)); %获取最接近的两个点的x坐标
            mode2nodes = mode(index(1:2), 1:mode_number); %获取两个点坐标的y值
            phi_result = interp1(xResult, mode2nodes, xTMD, 'linear', 'extrap'); %插值以后任意点的振型
            modeTMD(t1, 1:mode_number) = phi_result(1:mode_number);
        end

        %
        %     figure
        %     hold on
        %
        %     for k1 = 1:calmodes
        %         plot(nodegap, mode_re(:, k1))
        %     end
        %     legend('1', '2', '3')

        %     [modemaxdis_single,usinglemax,uallmax]= CalData_Polynomial_withTMD_multidegree(nTMD,mTMD,zetaTMD,omegaTMD,xTMD,1,ifcalmode,MM_eq,KK_eq,calmodes,eig_val,eig_vec,t_length);

%         [result]=CalDamping_Polynomial_withTMD_multidegree(nTMD,mTMD,zetaTMD,omegaTMD,xTMD,1,ifcalmode,MM_eq,KK_eq,calmodes,eig_val,eig_vec);
        [result]=CalDamping_Polynomial_withTMD_multidegree_no_aerodamping(nTMD,mTMD,zetaTMD,omegaTMD,xTMD,1,ifcalmode,MM_eq,KK_eq,calmodes,eig_val,eig_vec);
        

end

