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
