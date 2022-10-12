%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%Date: 2022-10-09 13:10:29
%LastEditors: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%LastEditTime: 2022-10-12 15:08:50
%FilePath: \NonlinearScanlan\20221012考虑一阶模态与多阶区别\CalData_Polynomial_withTMD_multidegree3.m
%Description:考虑一阶模态与多阶区别
%
%Copyright (c) 2022 by xushengyichn 54436848+xushengyichn@users.noreply.github.com, All Rights Reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 求解安装TMD后的响应

clc; clear; close all % 清除记录
addpath("../函数/")
numberofTMD = 1; % 所需要优化的TMD的数量
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
matrix=load('matrix.mat');
K=matrix.K;
M=matrix.M;

calmodes = 1; %考虑模态数 Consider the number of modes
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
    MM_eq(k1,k1)=P_eq(k1, mode_vec, M); %P=parameters
end

for k1 = 1:calmodes
    KK_eq(k1,k1)=P_eq(k1, mode_vec, K); %P=parameters
end

% save data
% load data
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

end

for k1= 1:nModes
    mode_re(:,k1)=mode(:,k1)./max(abs(mode(:,k1)));
end
[phideckmax,location] = max(mode); %计算桥面每个模态的最大振型值
nodegap = importdata('nodegap.txt');
modes=[nodegap nodeondeck mode mode_re];
figure
hold on
for k1= 1:calmodes
    plot(nodegap,mode_re(:,k1))
end
legend('1','2','3')
%% 计算不安装TMD情况下各阶模态各点最大位移
nTMD = 0;
mTMD = 0;
zetaTMD = 0;
omegaTMD = 0;
nodeTMD = 0;
mode_numbers = 1:1:1;
ifcalmode = 3;

for mode_number = 1:length(mode_numbers)
    [modemaxdis_single_noTMD(mode_number), usinglemax_noTMD(:, mode_number), uallmax_noTMD(:, mode_number)] = CalData_Polynomial_withTMD_multidegree(nTMD, mTMD, zetaTMD, omegaTMD, nodeTMD, mode_number, ifcalmode, MM_eq, KK_eq, calmodes, eig_val, eig_vec);
end

clear nTMD mTMD zetaTMD omegaTMD nodeTMD mode_numbers ifcalmode
%%
% 设置TMD参数

mass_six_span = 10007779.7;
mu = 0.001/100;
mTMDall = mass_six_span * mu;

nTMD = numberofTMD;
mTMD=mTMDall/nTMD*ones(nTMD,1);
Ftmd = 0.8339;
zetaTMD = 0.06;
omegaTMD = 2 * pi * Ftmd;
% cTMD = [2 * mTMD(1) * omegatmds * 0.05];
% disp(cTMD)
% 计算桥面节点的振型向量
KMmapping = importmappingmatrix('KMatrix.mapping');
UYNode = sort(KMmapping.Node(KMmapping.DOF == 'UY'));

nodeondeck = importdata('nodeondeck.txt');
nodegap = importdata('nodegap.txt');
nodeondecknew = [];
intalledlocation = [];

for k1 = 1:length(nodeondeck)

    if mod(k1, 5) == 0
        nodeondecknew = [nodeondecknew nodeondeck(k1)];
        intalledlocation = [intalledlocation nodegap(k1)];
    end

end

mode_number=1;%仅考虑一阶气动力的情况
% mode_numbers = 1:1:1;
ifcalmode = 3;
nodeTMD=[1036];
[modemaxdis_single,usinglemax,uallmax,u1]=CalData_Polynomial_withTMD_multidegree(nTMD,mTMD,zetaTMD,omegaTMD,nodeTMD,mode_number,ifcalmode,MM_eq,KK_eq,calmodes,eig_val,eig_vec);


phi1036_1=modes(find(modes(:,2)==1036),3);
dis1036_1mode=max(abs(u1(1,:)))*phi1036_1;

clearvars -except dis1036_1mode nTMD mTMD zetaTMD omegaTMD nodeTMD mode_number ifcalmode  u1 phi1036_1 mode_number
%% 考虑两阶模态振动
numberofTMD = 1; % 所需要优化的TMD的数量
% 导入质量刚度矩阵
matrix=load('matrix.mat');
K=matrix.K;
M=matrix.M;

calmodes = 2; %考虑模态数 Consider the number of modes
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
    MM_eq(k1,k1)=P_eq(k1, mode_vec, M); %P=parameters
end

for k1 = 1:calmodes
    KK_eq(k1,k1)=P_eq(k1, mode_vec, K); %P=parameters
end

% save data
% load data
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

end

for k1= 1:nModes
    mode_re(:,k1)=mode(:,k1)./max(abs(mode(:,k1)));
end
[phideckmax,location] = max(mode); %计算桥面每个模态的最大振型值
nodegap = importdata('nodegap.txt');
modes=[nodegap nodeondeck mode mode_re];
figure
hold on
for k1= 1:calmodes
    plot(nodegap,mode_re(:,k1))
end
legend('1','2','3')


% mode_numbers = 1:1:2;
ifcalmode = 3;
nodeTMD=[1036];
[modemaxdis_single,usinglemax,uallmax,u2]=CalData_Polynomial_withTMD_multidegree(nTMD,mTMD,zetaTMD,omegaTMD,nodeTMD,mode_number,ifcalmode,MM_eq,KK_eq,calmodes,eig_val,eig_vec);

phi1036_1=modes(find(modes(:,2)==1036),3);
phi1036_2=modes(find(modes(:,2)==1036),4);
dis1036_2mode=max(abs(u2(1,:)*phi1036_1+u2(2,:)*phi1036_2));
h = 0.01; % Time step
t = 0:h:100; % Time
close all
figure

plot(t,u1(1,:)*phi1036_1)
hold on
plot(t,u2(1,:)*phi1036_1+u2(2,:)*phi1036_2)
data=u2(1,:)*phi1036_1+u2(2,:)*phi1036_2;

disp(dis1036_1mode)
disp(dis1036_2mode)
disp(max(u1(end,:)))
disp(max(u2(end,:)))

 [psd_avg, f, psd_plot] = fft_transfer(1/h,data');
figure
plot(f,psd_plot)
% 所需函数

function result = P_eq(mode, temp_vec, Matrix)
    vec = temp_vec(:, mode);
    result = vec' * Matrix * vec;
end
