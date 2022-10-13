%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%Date: 2022-10-13 15:27:10
%LastEditors: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%LastEditTime: 2022-10-13 16:42:50
%FilePath: \NonlinearScanlan\20221012考虑一阶模态与多阶区别\DrawMaterial.m
%Description: 画两阶模态示意图
%
%Copyright (c) 2022 by xushengyichn 54436848+xushengyichn@users.noreply.github.com, All Rights Reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;

% 导入质量刚度矩阵
matrix=load('matrix.mat');
K=matrix.K;
M=matrix.M;

calmodes = 5; %考虑模态数 Consider the number of modes
[eig_vec, eig_val] = eigs(K, M, calmodes, 'SM');


% 计算桥面节点的振型向量
nModes=calmodes;
KMmapping = importmappingmatrix('KMatrix.mapping');
nodeondeck = importdata('nodeondeck.txt');
mode = zeros(length(nodeondeck), nModes);

[omeg, w_order] = sort(sqrt(diag(eig_val)));
mode_vec = eig_vec(:, w_order);
Freq = omeg / (2 * pi);

for k1 = 1:length(nodeondeck)
    position_index = KMmapping.MatrixEqn(find(and(KMmapping.Node == nodeondeck(k1), KMmapping.DOF == 'UY')));

    if isempty(position_index)
        mode(k1, :) = zeros(1, nModes);
    else
        mode(k1, :) = mode_vec(position_index, :);
    end

end

[phideckmax,location] = max(mode); %计算桥面每个模态的最大振型值
nodegap = importdata('nodegap.txt');

for k1 = 1:calmodes
    mode_re(:,k1)=mode(:,k1)/max(abs(mode(:,k1)));
end

figure
plot(nodegap,mode_re(:,1))

modedata=[nodeondeck nodegap mode mode_re];
modeinfo=array2table(modedata,"VariableNames",["node_number","node_x","mode1","mode2","mode3","mode4","mode5","mode_re1","mode_re2","mode_re3","mode_re4","mode5_re"]);

save modedata modedata