%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%Date: 2022-10-15 21:56:39
%LastEditors: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%LastEditTime: 2022-10-15 22:11:05
%FilePath: \NonlinearScanlan\20221012考虑一阶模态与多阶区别\Cal_Dis.m
%Description: 计算某个工况的响应
%
%Copyright (c) 2022 by xushengyichn 54436848+xushengyichn@users.noreply.github.com, All Rights Reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 求解安装TMD后的响应

clc; clear; close all % 清除记录
addpath("函数/")

numberofTMD = 2; % 所需要优化的TMD的数量.

savedata=1;

nTMD = numberofTMD;
% nTMD = 2; %TMD数量

mass_six_span = 10007779.7;
mu = [0.015 0.015];
mTMD = mu*mass_six_span;
Ftmd = [0.82 0.93];
zetaTMD = [0.1 0.22];
xTMD = [56.2 289.5];

temp=importdata("opt_1TMD_1_2modes_xloc56.2.mat");
a_mTMD=temp.sol.mTMD;
a_xTMD=temp.sol.xTMD;
a_Ftmd=temp.sol.fTMD;
a_zetaTMD=temp.sol.zetaTMD;

temp=importdata("opt_1TMD_1modes_xloc56.2.mat");
b_mTMD=temp.sol.mTMD;
b_xTMD=temp.sol.xTMD;
b_Ftmd=temp.sol.fTMD;
b_zetaTMD=temp.sol.zetaTMD;

mTMD=[a_mTMD b_mTMD];
xTMD=[a_xTMD b_xTMD];
Ftmd=[a_Ftmd b_Ftmd];
zetaTMD=[a_zetaTMD b_zetaTMD];


omegaTMD = 2 * pi * Ftmd;



mode_numbers = 1:1:2;
ifcalmode = 3;
h=0.01;
t_length=500;

%%
% 读入ANSYS质量刚度矩阵

% 读取梁桥模型

% 将ANSYS中的稀疏矩阵处理为完全矩阵
% Handling sparse matrices in ANSYS as full matrices 导入ANSYS MCK矩阵
% Import MCK matrix from ANSYS
hb_to_mm ('KMatrix.matrix', 'K.txt');
hb_to_mm ('MMatrix.matrix', 'M.txt');
% hb_to_mm ('CMatrix.matrix', 'C.txt');

%map the node and matrix from the KMatrix.mapping and MMatrix.mapping
Kdata = importdata('K.txt').data;
Kmatrix = zeros(Kdata(1, 1), Kdata(1, 2));

for i = 2:size(Kdata, 1)
    Kmatrix(Kdata(i, 1), Kdata(i, 2)) = Kdata(i, 3);
end

Mdata = importdata('M.txt').data;
Mmatrix = zeros(Mdata(1, 1), Mdata(1, 2));

for i = 2:size(Mdata, 1)
    Mmatrix(Mdata(i, 1), Mdata(i, 2)) = Mdata(i, 3);
end

% Cdata = importdata('C.txt').data;
% Cmatrix_DP = zeros(Cdata(1, 1), Cdata(1, 2));
%
% for i = 2:size(Cdata, 1)
%     Cmatrix_DP(Cdata(i, 1), Cdata(i, 2)) = Cdata(i, 3);
% end

% 还原对角线以上元素，使之为对称阵, ANSYS只给出下三角矩阵
% Restore the elements above the diagonal to make it a symmetric matrix, ANSYS only gives the lower triangular matrix
K = diag(diag(Kmatrix) / 2) + Kmatrix - diag(diag(Kmatrix));
K = K + K';
M = diag(diag(Mmatrix) / 2) + Mmatrix - diag(diag(Mmatrix));
M = M + M';
% C_exp = diag(diag(Cmatrix_DP) / 2) + Cmatrix_DP - diag(diag(Cmatrix_DP));
% C_exp = C_exp + C_exp';
% C =C_exp;
%%
% 模态矩阵计算

% 导入质量刚度矩阵
matrix = load('matrix.mat');
K = matrix.K;
M = matrix.M;

%% 计算不安装TMD情况下各阶模态各点最大位移

t = 0:h:t_length; % Time
for mode_number = 1:length(mode_numbers)
    calmodes = mode_number; %考虑模态数 Consider the number of modes
    [eig_vec, eig_val] = eigs(K, M, calmodes, 'SM');

    [nfdof, nfdof] = size(eig_vec);

    for j = 1:nfdof
        mnorm = sqrt(eig_vec(:, j)' * M * eig_vec(:, j));
        eig_vec(:, j) = eig_vec(:, j) / mnorm; %振型质量归一化 Mode shape mass normalization
    end

    for j = 1:nfdof
        m_modal(j) = sqrt(eig_vec(:, j)' * M * eig_vec(:, j));
    end

    [omeg, w_order] = sort(sqrt(diag(eig_val)));
    mode_vec = eig_vec(:, w_order);
    Freq = omeg / (2 * pi);

    for k1 = 1:calmodes
        MM_eq(k1, k1) = P_eq(k1, mode_vec, M); %P=parameters
    end

    for k1 = 1:calmodes
        KK_eq(k1, k1) = P_eq(k1, mode_vec, K); %P=parameters
    end

    % save data
    % load data
    nModes = calmodes;
    % 计算桥面节点的振型向量
    nodeondeck = importdata('nodeondeck.txt');
    KMmapping = importmappingmatrix('KMatrix.mapping');
    UYNode = sort(KMmapping.Node(KMmapping.DOF == 'UY'));
    mode = zeros(length(nodeondeck), nModes);

    for k1 = 1:length(nodeondeck)
        position_index = KMmapping.MatrixEqn(find(and(KMmapping.Node == nodeondeck(k1), KMmapping.DOF == 'UY')));

        if isempty(position_index)
            mode(k1, :) = zeros(1, nModes);
        else
            mode(k1, :) = mode_vec(position_index, :);
        end

    end

    for k1 = 1:nModes
        mode_re(:, k1) = mode(:, k1) ./ max(abs(mode(:, k1)));
    end

    [phideckmax, location] = max(mode); %计算桥面每个模态的最大振型值
    nodegap = importdata('nodegap.txt');
    modes = [nodegap nodeondeck mode mode_re];
%      if length(mTMD)==1
%     for k1 = 1:length(mTMD)
%         modeTMD(k1,1:mode_number)=mode(nodeTMD(k1)==modes(:,2),:);
%     end

for t1 = 1:length(mTMD)

        [~,index]=sort(abs(nodegap-xTMD));%查找与xTMD最接近的点的排序
        xResult=nodegap(index(1:2));%获取最接近的两个点的x坐标
        mode2nodes=mode(index(1:2),1:mode_number);%获取两个点坐标的y值
        phi_result=interp1(xResult,mode2nodes,xTMD,'linear','extrap');%插值以后任意点的振型
        modeTMD(t1, 1:mode_number) = phi_result(1:mode_number);
end

%      end
    figure
    hold on

    for k1 = 1:calmodes
        plot(nodegap, mode_re(:, k1))
    end

    legend('1', '2', '3')
%     [modemaxdis_single,usinglemax,uallmax]= CalData_Polynomial_withTMD_multidegree(nTMD,mTMD,zetaTMD,omegaTMD,xTMD,1,ifcalmode,MM_eq,KK_eq,calmodes,eig_val,eig_vec,t_length);
  [modemaxdis_single,usinglemax,uallmax,output]=CalData_Polynomial_withTMD_multidegree(nTMD,mTMD,zetaTMD,omegaTMD,xTMD,1,ifcalmode,MM_eq,KK_eq,calmodes,eig_val,eig_vec);
  u=output.u;   
  if length(mTMD)==1%安装1个TMD的情况
    if mode_number==1
        uNodeTMD(mode_number,:)=modeTMD'.*u(1:end-nTMD,:);
    else
        uNodeTMD(mode_number,:)=sum(modeTMD'.*u(1:end-nTMD,:));
    end
    end

    if length(mTMD)==2%安装2个TMD的情况
        if mode_number==1
        uNodeTMD(mode_number,:)=modeTMD(2,1)*u(1:1,:);
        end
        if mode_number==2
        uNodeTMD(mode_number,:)=modeTMD(2,1)*u(1,:)+modeTMD(2,2)*u(2,:);
        end
    end
end
   if length(mTMD)==1
figure
plot(t,uNodeTMD(1,:))
hold on
plot(t,uNodeTMD(2,:))
   end
% fs=1/h;
% [psd_avg, f, psd_plot] = fft_transfer(fs,uNodeTMD(1,:)');
% figure
% plot(f, psd_plot)
% 
% [psd_avg, f, psd_plot] = fft_transfer(fs,uNodeTMD(2,:)');
% figure
% plot(f, psd_plot)
if savedata==1
    save uNode56.2_290With2TMD.mat uNodeTMD  %储存TMD安装在1036节点时的响应结果（分别考虑一阶或两阶模态的结果）
end
%% 所需函数

function result = P_eq(mode, temp_vec, Matrix)
    vec = temp_vec(:, mode);
    result = vec' * Matrix * vec;
end
