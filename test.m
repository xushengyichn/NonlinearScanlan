%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%Date: 2022-11-04 09:00:33
%LastEditors: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%LastEditTime: 2022-11-04 09:01:00
%FilePath: \NonlinearScanlan\test.m
%Description: 
%
%Copyright (c) 2022 by xushengyichn 54436848+xushengyichn@users.noreply.github.com, All Rights Reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear
close all

addpath("函数\")

nodegap = importdata('nodegap.txt');
nodeondeck = importdata('nodeondeck.txt');

load modeinfo.mat
nModes=5;
mode = zeros(length(nodeondeck), nModes);
KMmapping = importmappingmatrix('KMatrix.mapping');
for k1 = 1:length(nodeondeck)
    position_index = KMmapping.MatrixEqn(find(and(KMmapping.Node == nodeondeck(k1), KMmapping.DOF == 'UY')));

    if isempty(position_index)
        mode(k1, :) = zeros(1, nModes);
    else
        mode(k1, :) = mode_vec(position_index, :);
    end

end

xTMD=289.465;
% xTMD=56.2;
[~,index]=sort(abs(nodegap-xTMD));%查找与xTMD最接近的点的排序
xResult=nodegap(index(1:2));%获取最接近的两个点的x坐标
mode2nodes=mode(index(1:2),1:nModes);%获取两个点坐标的y值
phi_result=interp1(xResult,mode2nodes,xTMD);%插值以后任意点的振型


mb=1;
xib=0.05;
fb=1;
wb=2*pi*fb;
mu=0.01;
tau=1.05;
miu=10;

M=mb/mu*[1 0 ;0 mu];
C=2*wb*xib*mb/mu*[1/tau/miu+mu -mu;-mu mu];
K=wb^2*mb/mu*[1/tau/tau+mu -mu ;-mu mu];

disp(M)
disp(C)
disp(K)

Mode = Complex_Eigenvalue_Analysis(M,K*0.0001,K)

