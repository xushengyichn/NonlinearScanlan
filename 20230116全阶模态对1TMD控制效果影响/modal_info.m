%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%Date: 2022-10-17 13:50:41
%LastEditors: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%LastEditTime: 2022-10-17 13:53:07
%FilePath: \NonlinearScanlan\modal_info.m
%Description: 提取连续梁模态信息
%
%Copyright (c) 2022 by xushengyichn 54436848+xushengyichn@users.noreply.github.com, All Rights Reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc;clear;close all;
tic
% % 读入ANSYS质量刚度矩阵
% 
% % 读取梁桥模型
% 
% % 将ANSYS中的稀疏矩阵处理为完全矩阵
% % Handling sparse matrices in ANSYS as full matrices 导入ANSYS MCK矩阵
% % Import MCK matrix from ANSYS
% hb_to_mm ('KMatrix.matrix', 'K.txt');
% hb_to_mm ('MMatrix.matrix', 'M.txt');
% % hb_to_mm ('CMatrix.matrix', 'C.txt');
% 
% %map the node and matrix from the KMatrix.mapping and MMatrix.mapping
% Kdata = importdata('K.txt').data;
% Kmatrix = zeros(Kdata(1, 1), Kdata(1, 2));
% 
% for i = 2:size(Kdata, 1)
%     Kmatrix(Kdata(i, 1), Kdata(i, 2)) = Kdata(i, 3);
% end
% 
% Mdata = importdata('M.txt').data;
% Mmatrix = zeros(Mdata(1, 1), Mdata(1, 2));
% 
% for i = 2:size(Mdata, 1)
%     Mmatrix(Mdata(i, 1), Mdata(i, 2)) = Mdata(i, 3);
% end
% 
% % Cdata = importdata('C.txt').data;
% % Cmatrix_DP = zeros(Cdata(1, 1), Cdata(1, 2));
% %
% % for i = 2:size(Cdata, 1)
% %     Cmatrix_DP(Cdata(i, 1), Cdata(i, 2)) = Cdata(i, 3);
% % end
% 
% % 还原对角线以上元素，使之为对称阵, ANSYS只给出下三角矩阵
% % Restore the elements above the diagonal to make it a symmetric matrix, ANSYS only gives the lower triangular matrix
% K = diag(diag(Kmatrix) / 2) + Kmatrix - diag(diag(Kmatrix));
% K = K + K';
% M = diag(diag(Mmatrix) / 2) + Mmatrix - diag(diag(Mmatrix));
% M = M + M';
% % C_exp = diag(diag(Cmatrix_DP) / 2) + Cmatrix_DP - diag(diag(Cmatrix_DP));
% % C_exp = C_exp + C_exp';
% % C =C_exp;

addpath("../函数/")
% 导入矩阵

matrix=load('matrix.mat');
K=matrix.K;
M=matrix.M;
% M_inv=inv(M);
calmodes=100; % 计算模态数
[eig_vec, eig_val] = eigs(K, M, calmodes, 'SM');
% [eig_vec, eig_val] = eigs(K, M);

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
    MM_eq(k1,k1)=P_eq(k1, mode_vec, M); %P=parameters %模态质量矩阵
end

for k1 = 1:calmodes
    KK_eq(k1,k1)=P_eq(k1, mode_vec, K); %P=parameters %模态刚度矩阵
end



nModes=calmodes;
% 计算桥面节点的振型向量
nodeondeck = importdata('nodeondeck.txt');
KMmapping = importmappingmatrix('KMatrix.mapping');
UYNode=sort(KMmapping.Node(KMmapping.DOF == 'UY'));
mode = zeros(length(nodeondeck), nModes);

for k1 = 1:length(nodeondeck)
    position_index = KMmapping.MatrixEqn(find(and(KMmapping.Node == nodeondeck(k1), KMmapping.DOF == 'UY')));

    if isempty(position_index)
        mode(k1, :) = zeros(1, nModes);
    else
        mode(k1, :) = mode_vec(position_index, :);
    end
    clear position_index
end

for k1= 1:nModes
    mode_re(:,k1)=mode(:,k1)./max(abs(mode(:,k1)));
end
[phideckmax,location] = max(mode); %计算桥面每个模态的最大振型值
nodegap = importdata('nodegap.txt');
modes=[nodegap nodeondeck mode mode_re];

clear K M calmodes j k1 matrix mnorm nfdof nModes 
save modeinfo_all 
toc
function result = P_eq(mode, temp_vec, Matrix)
    vec = temp_vec(:, mode);
    result = vec' * Matrix * vec;
end