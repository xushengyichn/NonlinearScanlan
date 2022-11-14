%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%Date: 2022-10-15 21:56:39
%LastEditors: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%LastEditTime: 2022-11-14 23:07:46
%FilePath: \NonlinearScanlan\20221114一阶模态二TMD结果穷举\Cal_Dis_onemode_twoTMD_zeta_fre_loc_three_parameters.m
%Description: 计算一阶模态，两个TMD的影响，穷举阻尼频率和位置
%
%Copyright (c) 2022 by xushengyichn 54436848+xushengyichn@users.noreply.github.com, All Rights Reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 求解安装TMD后的响应

clc; clear; close all % 清除记录
addpath("../函数/")

numberofTMD = 2; % 所需要计算的TMD的数量.

savedata = 0;

nTMD = numberofTMD;
% nTMD = 2; %TMD数量
modeinfo = load('modeinfo.mat');
% 设计第一个TMD的参数
% TMD1按照规范设计
fs=modeinfo.Freq(1);
mode1 = modeinfo.eig_vec(:, 1);
mu = 0.02;
mTMD1 = mu / (max(mode1)^2);
fTMD1 = 1/(1+mu)*fs;
zetaTMD1 = sqrt(3*mu/8/(1+mu));

xTMD2_all=0:1:660;
fTMD2_all=fTMD1*0.8:0.01:fTMD1*1.2;
zetaTMD2_all=zetaTMD1*0.5:0.01:zetaTMD1*3;

[XTMD2_all,FTMD2_all]=ndgrid(xTMD2_all,fTMD2_all);
variables = [XTMD2_all(:),FTMD2_all(:)];

onemode_twotmd_dis=zeros(size(variables,1),4);%第四列为是否收敛标志
numIterations=size(variables,1);

% ppm = ParforProgressbar(numIterations,'showWorkerProgress',true,'progressBarUpdatePeriod',3,'title','my fancy title'); 


pauseTime = 60/numIterations;
% parfor k1 = 1:size(variables,1)
for k1 = 1:5

mass_six_span = 10007779.7;
mTMD = [mTMD1 mTMD1];
Ftmd = [fTMD1 variables(k1,2)];
zetaTMD = [zetaTMD1 zetaTMD1];
xTMD = [384.5 variables(k1,1)];

omegaTMD = 2 * pi * Ftmd;

mode_numbers = 1;
ifcalmode = 3;
h = 0.01;
t_length =150;

%% 计算不安装TMD情况下各阶模态各点最大位移
nodeondeck = importdata('nodeondeck.txt');
KMmapping = importmappingmatrix('KMatrix.mapping');
t = 0:h:t_length; % Time

for mode_number = 1:length(mode_numbers)
    iter=1;
    flag_convergence=0;
    while and(iter<=5,flag_convergence==0)
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
    
    mode_re=zeros(size(mode,1),size(mode,2));
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
    
    modeTMD=zeros(length(mTMD),mode_number);
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
    [modemaxdis_single, usinglemax, uallmax, output] = CalData_Polynomial_withTMD_multidegree(nTMD, mTMD, zetaTMD, omegaTMD, xTMD, 1, ifcalmode, MM_eq, KK_eq, calmodes, eig_val, eig_vec,t_length);
    u = output.u;

    dis_beam=max(mode1)*u(1,:);
    dis_TMD1=u(2,:);
    dis_TMD2=u(3,:);

    %判断计算是否收敛
    dis1=max(dis_beam(end/8*6:end/8*7));
    dis2=max(dis_beam(end/8*7:end));
    
        if or(abs((dis1-dis2)/dis2)<0.05,dis2<1e-3)
            disp("计算收敛")
            flag_convergence=1;
            dis_beam_max=std(dis_beam(end/8*6:end))*sqrt(2);
            dis_TMD1=std(u(2,end/8*6:end))*sqrt(2);
            dis_TMD2=std(u(3,end/8*6:end))*sqrt(2); 
            tnew=0:0.01:t_length;
            if iter>=3
                flag_iter=0;
            else
                flag_iter=1;
            end
            result_temp=[dis_beam_max dis_TMD1 dis_TMD2 flag_iter];
            onemode_twotmd_dis(k1,:)=result_temp;
%             figure
%             plot(tnew,max(mode1)*u(1,:))   
        else
            disp("计算未收敛，增加计算时间")
            iter=iter+1;
            tnew=0:0.01:t_length;
%             figure
%             plot(tnew,max(mode1)*u(1,:))   
            t_length=t_length+300;
        end
    end
    

    
end
%     pause(pauseTime);
%     % increment counter to track progress
%     ppm.increment();
end

onemode_twotmd_results=[variables onemode_twotmd_dis];


[XTMD_all,FTMD2_all]=ndgrid(xTMD2_all,fTMD2_all);
variables = [XTMD_all(:),FTMD2_all(:)];
bridge_dis_grid=griddata(variables(:,1),variables(:,2),onemode_twotmd_dis(:,1),xTMD2_all,FTMD2_all);
TMD1_dis_grid=griddata(variables(:,1),variables(:,2),onemode_twotmd_dis(:,2),xTMD2_all,FTMD2_all);
TMD2_dis_grid=griddata(variables(:,1),variables(:,2),onemode_twotmd_dis(:,3),xTMD2_all,FTMD2_all);

save onemode_twotmd_results.mat onemode_twotmd_phi_results xTMD2_all FTMD2_all bridge_dis_grid TMD1_dis_grid TMD2_dis_grid

% save onemode_twotmd_results.mat onemode_twotmd_results

%% 所需函数

function result = P_eq(mode, temp_vec, Matrix)
    vec = temp_vec(:, mode);
    result = vec' * Matrix * vec;
end
