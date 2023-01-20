%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%Date: 2022-10-15 21:56:39
%LastEditors: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%LastEditTime: 2023-01-20 11:12:29
%FilePath: \NonlinearScanlan\20230116全阶模态对1TMD控制效果影响\Cal_Dis_allmodes_oneTMD_loc_two_parameters_aerofoce_17modes.m
%Description: 由于前8阶模态似乎不够，需要考虑更多阶模态，本代码修改为17阶模态（其中包含12阶竖弯模态)
%
%Copyright (c) 2022 by xushengyichn 54436848+xushengyichn@users.noreply.github.com, All Rights Reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 求解安装TMD后的响应

clc; clear; close all % 清除记录
addpath("../函数/")

airmodes = 1:5;

for a1 = airmodes %气动力作用在前5阶模态上

    numberofTMD = 1; % 所需要计算的TMD的数量.

    savedata = 0; % 是否保存数据

    nTMD = numberofTMD; %TMD数量
    % nTMD = 2; %TMD数量
    modeinfo = load('modeinfo_all.mat');
    
    modes_number = 17; % 求解前n阶模态
    fs = modeinfo.Freq(a1); % 按照气动力作用在哪个模态上，来确定TMD的频率

    % 获取前n阶模态的模态位移，这里的mode_shape是按照桥梁模态逐阶排列的，不是按照cal_modes_index排列的，重新排列的向量为后续的eig_vec变量
    mode_shape = zeros(length(modeinfo.eig_vec(:, a1)), modes_number);

    for k1 = 1:modes_number
        mode_shape(:, k1) = modeinfo.eig_vec(:, k1);
    end

    % 设计第一个TMD的参数
    % TMD1按照规范设计，按照气动力对应模态的最优参数进行设计
    mu = 0.02;
    mTMD1 = mu / (max(mode_shape(:, a1)) ^ 2);
    fTMD1 = 1 / (1 + mu) * fs;
    zetaTMD1 = sqrt(3 * mu / 8 / (1 + mu));

    cal_modes = (1:modes_number)';% 

    % cal_modes = [cal_modes; 100];

    % cal_modes_index 先计算对应气动力的单模态响应,再由小到大逐步增加模态数直至达到对应要求的模态数modes_number.代码逻辑为根据cal_modes的值,从cal_modes_index来查找考虑的实际模态
    cal_modes_index=(1:modes_number)';
    temp= find(cal_modes_index==a1);
    cal_modes_index(temp)=[];
    cal_modes_index=[a1;cal_modes_index];

    cal_modes_accurate=100; % 精确计算的模态数
    % TMD1的位置
    xTMD1_all = 0:1:660;

    [XTMD1_all, Cal_modes] = ndgrid(xTMD1_all, cal_modes);
    variables = [XTMD1_all(:), Cal_modes(:)];

    nmodes_onetmd_dis = zeros(size(variables, 1), 4); %四列分别代表，最后一列为是否收敛dis_beam_max dis_TMD1_max max_index flag_iter
    numIterations = size(variables, 1);

    ppm = ParforProgressbar(numIterations, 'showWorkerProgress', true, 'progressBarUpdatePeriod', 3, 'title', 'my fancy title');

    pauseTime = 60 / numIterations;

    parfor k1 = 1:size(variables, 1)
        % parfor k1 = 1:modes_number
        % for k1 = 537:603
        % for k1 = 63
    % for k1 = 1267

        % mass_six_span = 10007779.7;
        mTMD = [mTMD1]; %该代码为基准tmd的
        Ftmd = [fTMD1];
        zetaTMD = [zetaTMD1];
        xTMD = [variables(k1, 1)];

        omegaTMD = 2 * pi * Ftmd;

        cal_modes_temp = variables(k1, 2);
        cal_modes_index_temp=cal_modes_index(1:cal_modes_temp);


        ifcalmode = 3;
        h = 0.01;
        t_length = 300;

        %% 计算不安装TMD情况下各阶模态各点最大位移
        nodeondeck = importdata('nodeondeck.txt');
        KMmapping = importmappingmatrix('KMatrix.mapping');
        t = 0:h:t_length; % Time

%         mode_number = mode_numbers;
        iter = 1;
        flag_convergence = 0;

        while and(iter <= 3, flag_convergence == 0)
            calmodes = cal_modes_temp; %考虑模态数 Consider the number of modes
            nModes = calmodes;
            
            % 生成计算用质量刚度矩阵以及模态信息
            MM_eq = zeros(calmodes, calmodes);
            KK_eq = zeros(calmodes, calmodes);
            eig_val = zeros(calmodes, calmodes);
            eig_vec = zeros(length(modeinfo.eig_vec(:, 1)), calmodes);
            % 将气动力施加的模态始终放在第一位，随后补充其余模态
            for k2 = 1:length(cal_modes_index_temp)
                MM_eq(k2,k2)=modeinfo.MM_eq(cal_modes_index_temp(k2),cal_modes_index_temp(k2));
                KK_eq(k2,k2)=modeinfo.KK_eq(cal_modes_index_temp(k2),cal_modes_index_temp(k2));
                eig_val(k2,k2)=modeinfo.eig_val(cal_modes_index_temp(k2),cal_modes_index_temp(k2));
                eig_vec(:,k2)=modeinfo.eig_vec(:,cal_modes_index_temp(k2));
            end
            % MM_eq = modeinfo.MM_eq(1:calmodes, 1:calmodes);
            % KK_eq = modeinfo.KK_eq(1:calmodes, 1:calmodes);
            % eig_val = modeinfo.eig_val(1:calmodes, 1:calmodes);
            % eig_vec = modeinfo.eig_vec(:, 1:calmodes);

            % 计算桥面节点的振型向量

            UYNode = sort(KMmapping.Node(KMmapping.DOF == 'UY'));
            mode = zeros(length(nodeondeck), nModes); %计算桥面的模态数
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
            %         modeTMD(k1,1:mode_number)=mode(nodeTMD(k1)==modes(:,2),:);dd
            %     end

            modeTMD = zeros(length(mTMD), calmodes);

            for t1 = 1:length(mTMD)

                [~, index] = sort(abs(nodegap - xTMD)); %查找与xTMD最接近的点的排序
                xResult = nodegap(index(1:2)); %获取最接近的两个点的x坐标
                mode2nodes = mode(index(1:2), 1:calmodes); %获取两个点坐标的y值
                phi_result = interp1(xResult, mode2nodes, xTMD, 'linear', 'extrap'); %插值以后任意点的振型
                modeTMD(t1, 1:calmodes) = phi_result(1:calmodes);
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
            [~, ~, ~, output] = CalData_Polynomial_withTMD_multidegree_8modes(nTMD, mTMD, zetaTMD, omegaTMD, xTMD,1, ifcalmode, MM_eq, KK_eq, calmodes, eig_val, eig_vec, t_length);
            u = output.u;

            %计算桥梁响应
            dis_temp = zeros(size(u, 2), variables(k1, 2));% 分别表示时间长度和模态数量

            for k4 = 1:variables(k1, 2)
                dis_temp(:, k4) = u(k4, :);
            end

            maxdislocationtime = u(1, :);

            dis_temp(:, variables(k1, 2) + 1) = u(end, :); %TMD响应
            dis_TMD = dis_temp(:, variables(k1, 2) + 1);
            %判断计算是否收敛
            dis1 = max(maxdislocationtime(round(end / 8 * 6, 0):round(end / 8 * 7, 0)));
            dis2 = max(maxdislocationtime(round(end / 8 * 7, 0):end));

            if or(abs((dis1 - dis2) / dis2) < 0.02, dis2 < 1e-3)
                disp("计算收敛")
                flag_convergence = 1;

                npoints = size(eig_vec, 1);
                seq_length = length(u(1, :));
                nModes = nModes;

                %             seqs_all=zeros(npoints,seq_length,nModes);
                seqs_allmodes = zeros(npoints, seq_length);
                %             for k6=1:nModes
                %                 seqs_all(:,:,k6)=mode_shape(:,k6).*u(k6,:);
                %             end

                for k6 = 1:nModes
                    seqs_allmodes = seqs_allmodes + eig_vec(:, k6) .* u(k6, :);
                end

                %             seqs_1modes_max=max(seqs_all(:,:,1),[],2);
                %             seqs_2modes_max=max(seqs_all(:,:,2),[],2);
                %             seqs_3modes_max=max(seqs_all(:,:,3),[],2);
                %             seqs_4modes_max=max(seqs_all(:,:,4),[],2);
                %             seqs_5modes_max=max(seqs_all(:,:,5),[],2);

                seqs_allmodes_max = max(seqs_allmodes, [], 2);
                %寻找最大值和位置
                [max_value, max_index] = max(seqs_allmodes_max);
                seqs_allmodes_max_point = seqs_allmodes(max_index, :);
                seqs_allmodes = [];
                %             seqs_all=[];
                dis_beam_max = std(seqs_allmodes_max_point(round(end / 8 * 6, 0):end)) * sqrt(2);
                dis_TMD1_max = std(dis_TMD(round(end / 8 * 6, 0):end)) * sqrt(2);
                %             dis_TMD2=std(u(3,end/8*6:end))*sqrt(2);
                tnew = 0:0.01:t_length;

                if iter >= 3
                    flag_iter = 0;
                else
                    flag_iter = 1;
                end

                result_temp = [dis_beam_max dis_TMD1_max max_index flag_iter];
                nmodes_onetmd_dis(k1, :) = result_temp;
                %             figure
                %             plot(tnew,max(mode1)*u(1,:))
            else
                disp("计算未收敛，增加计算时间")
                iter = iter + 1;
                tnew = 0:0.01:t_length;
                %             figure
                %             plot(tnew,max(mode1)*u(1,:))
                t_length = t_length * 2 + 150;
            end

        end

%         pause(pauseTime);
%         % increment counter to track progress
%         ppm.increment();

    end

%     delete(ppm)
    nmodes_onetmd_results_loc = [variables nmodes_onetmd_dis];
    collectdata.nmodes_onetmd_results_loc=nmodes_onetmd_results_loc;
    collectdata.calmodes_index=cal_modes_index;
    str = "save('17modes_onetmd_results_loc_mode"+num2str(airmodes(a1))+".mat', 'collectdata')";
    eval(str);
    % save 10modes_onetmd_results_loc.mat nmodes_onetmd_results_loc
end %airmodes

%% 所需函数

function result = P_eq(mode, temp_vec, Matrix)
    vec = temp_vec(:, mode);
    result = vec' * Matrix * vec;
end
