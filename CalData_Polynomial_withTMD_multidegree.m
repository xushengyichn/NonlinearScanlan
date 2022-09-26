%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%Date: 2022-09-26 19:35:04
%LastEditors: Shengyi Xu xushengyichn@outlook.com
%LastEditTime: 2022-09-27 00:27:46
%FilePath: \NonlinearScanlan\CalData_Polynomial_withTMD_multidegree.m
%Description: 计算多模态，施加某一阶模态多项式气动力模型后的响应，考虑TMD
%
%Copyright (c) 2022 by xushengyichn 54436848+xushengyichn@users.noreply.github.com, All Rights Reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;

%% 参数设置

D=20; %断面参考宽度

nTMD = 1;
mTMD = [2400];
cTMD = [2 * mTMD(1) * 0.069704171453635 * 0.05];
kTMD = [mTMD(1) * 0.069704171453635];
nodeTMD = [2157]; %Node number(location of the TMD)

girderindex=1;

ExpNames = [
    'SZTD-110-case2-22.3-fasan-2401';
    'SZTD-110-case2-22.3-fasan-2501';
    'SZTD-110-case2-22.3-fasan-2601';
    'SZTD-110-case2-22.3-fasan-2701';
    'SZTD-110-case2-22.3-fasan-2801';
    'SZTD-110-case2-22.3-fasan-2901';
    'SZTD-110-case2-22.3-fasan-3101';
    ]; %记录文件名
    
% 设置模态叠加法生成阻尼矩阵的方式
% 1 表示采用瑞利阻尼矩阵生成对应的模态叠加法阻尼矩阵
% 2 表示采用模态阻尼比直接生成的模态叠加法阻尼矩阵，需要指定模态阻尼比zeta
% Set the method for generating the damping matrix by the modal superposition method
% 1 means use the Rayleigh damping matrix to generate the corresponding modal superposition method damping matrix
% 2 means use the modal damping ratio to directly generate the modal superposition method damping matrix, you need to specify the modal Damping ratio zeta



Zeta0=0.3/100;
disp("模态阻尼比未定义，采用默认值，每阶模态阻尼比均为："+num2str(Zeta0))

mode_number=1; %气动力施加的模态


%% 读取梁桥模型

%% 将ANSYS中的稀疏矩阵处理为完全矩阵
%% Handling sparse matrices in ANSYS as full matrices
% 导入ANSYS MCK矩阵
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


% 设置模态叠加法生成阻尼矩阵的方式
% 1 表示采用瑞利阻尼矩阵生成对应的模态叠加法阻尼矩阵
% 2 表示采用模态阻尼比直接生成的模态叠加法阻尼矩阵，需要指定模态阻尼比zeta
% Set the method for generating the damping matrix by the modal superposition method
% 1 means use the Rayleigh damping matrix to generate the corresponding modal superposition method damping matrix
% 2 means use the modal damping ratio to directly generate the modal superposition method damping matrix, you need to specify the modal Damping ratio zeta
DampingMatrixParameter = 2;

% 设置桥面坐标导入模式
% 1 采用程序自动写入一个含有桥面节点的数组，从 1 到 n
% 2 读入文件
% Set the bridge deck coordinate import mode
% 1 Use the program to automatically write an array containing the bridge deck nodes, from 1 to n
% 2 read in file
nodeondeckimport = 1;
% 设置外荷载施加方式
% 1 施加模态力
% 2 施加节点力（需要调整对应施加荷载部分的代码，未完成）
% Set the external load application method
% 1 Apply modal force
% 2 Apply nodal forces (need to adjust the code corresponding to the applied load section)
externalforcemethod = 1;
% 是否输出振动视频
% Whether to output vibration video
% 1 export vides
% 2 do not export vides
exportvideo=2;


% 特征值分析，即计算频率Freq和振型Phi，calmodes数字代表求解的阶数，eigs中参数SM表示从较小的特征值开始求解
% Eigenvalue analysis, that is to calculate the frequency Freq and mode shape Phi, the calmodes number represents the order of the solution, and the parameter SM in eigs represents the solution from the smaller eigenvalue.
calmodes = 10; %考虑模态数 Consider the number of modes
nModes=calmodes;
[eig_vec, eig_val] = eigs(K, M, calmodes, 'SM');
[nfdof, nfdof] = size(eig_vec);

for j = 1:nfdof
    mnorm = sqrt(eig_vec(:, j)' * M * eig_vec(:, j));
    eig_vec(:, j) = eig_vec(:, j) / mnorm; %振型质量归一化 Mode shape mass normalization
end

[omeg, w_order] = sort(sqrt(diag(eig_val)));
mode_vec = eig_vec(:, w_order);
Freq = omeg / (2 * pi);




%% 构建包含TMD的MCK矩阵
% 导入矩阵序号对应节点自由度关系。导入KCM中 *.mapping 文件的任意一个即可。他们是一样的。
% The imported matrix number corresponds to the node degree of freedom relationship. Import any one of the *.mapping files in KCM. They are the same.
KMmapping = importmappingmatrix('KMatrix.mapping');

% 计算桥面节点的振型向量
nodeondeck=importdata('nodeondeck.txt');
mode=zeros(length(nodeondeck),nModes);
for k1 = 1:length(nodeondeck)
    position_index = KMmapping.MatrixEqn(find(and(KMmapping.Node == nodeondeck(k1), KMmapping.DOF == 'UY')));
    if isempty(position_index)
        mode(k1,:) = zeros(1,nModes);
    else
        mode(k1,:) = mode_vec(position_index,:);
    end
end

nodegap=importdata('nodegap.txt');
for k1 =1:nModes
    integral_2(k1)=sum(mode(:,k1).^2.*nodegap);
    integral_3(k1)=sum(mode(:,k1).^2.*abs(mode(:,k1)).*nodegap);
    integral_4(k1)=sum(mode(:,k1).^4.*nodegap);
    integral_5(k1)=sum(mode(:,k1).^2.*abs(mode(:,k1)).^3.*nodegap);
    integral_6(k1)=sum(mode(:,k1).^6.*nodegap);
end


phi=mode(:,mode_number);
mode_integral_2=integral_2(mode_number);
mode_integral_3=integral_3(mode_number);
mode_integral_4=integral_4(mode_number);
mode_integral_5=integral_5(mode_number);
mode_integral_6=integral_6(mode_number);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 考虑TMD振动响应的模态叠加法
%% Modal Superposition Method Considering TMD Vibration Response

matrixsize = nTMD + nModes;
phiTMD = zeros(nTMD, nModes);
% phiTMD row:TMD for each loaction column:the mode shape at the each
% location of tmd
for t1 = 1:nTMD

    for t2 = 1:nModes
        % 找到TMD安装位置对应的振型向量大小并储存到phiTMD变量中
        % Find the size of the mode shape vector corresponding to the installation position of the TMD and store it in the phiTMD variable
        position_index = KMmapping.MatrixEqn(find(and(KMmapping.Node == nodeTMD(t1), KMmapping.DOF == 'UY')));
        phiTMD(t1, t2) = mode_vec(position_index, t2);
    end

end

clear t1 t2

% 创建质量矩阵
% Create the Mass matrix
MM = zeros(matrixsize, matrixsize);

for t1 = 1:matrixsize

    if t1 <= nModes
        MM(t1, t1) = P_eq(t1, mode_vec, M); %行列数小于等于nModes为模态质量 The number of rows and columns is less than nModes is the modal quality
    end

    if and(t1 > nModes, t1 <= matrixsize)
        MM(t1, t1) = mTMD(t1 - nModes); %行列数大于nModes为TMD实际质量 The number of rows and columns is greater than nModes for the actual quality of TMD
    end

end

clear t1

% 创建刚度矩阵
% Create the Stiffness Matrix
KK = zeros(size(MM, 1), size(MM, 2));
KK1 = zeros(size(MM, 1), size(MM, 2));
KK2 = zeros(size(MM, 1), size(MM, 2));
for k1 = 1:matrixsize
        if k1<=nModes
            KK1(k1,k1)=MM(k1,k1)*(2*pi*Freq(k1))^2;
        elseif k1>nModes
            KK1(k1,k1)=kTMD(k1-nModes);
        end
end
for k1 = 1:matrixsize
    if k1<=nModes
        for k2 = 1:length(mTMD)
            KK2(k1,k1)=KK2(k1,k1)+kTMD(k2);
        end
    elseif k1>nModes
        KK2(1,k1)=-kTMD(k1-nModes);
        KK2(k1,1)=-kTMD(k1-nModes);
    end
end
KK = KK1 + KK2;
clear k1 k2

% 创建阻尼矩阵
% Create the Damping Matrix
CC  = zeros(size(MM, 1), size(MM, 2));
CC1 = zeros(size(MM, 1), size(MM, 2));
CC2 = zeros(size(MM, 1), size(MM, 2));

for k1 = 1:matrixsize

        if k1<=nModes
            CC1(k1,k1)=Zeta0 * 4 * pi * MM(k1,k1) * Freq(k1);
        elseif k1>nModes
            CC1(k1,k1)=cTMD(k1-nModes);
        end
end

for k1 = 1:matrixsize
    if k1<=nModes
        for k2 = 1:length(mTMD)
            CC2(k1,k1)=CC2(k1,k1)+cTMD(k2);
        end
    elseif k1>nModes
        CC2(1,k1)=-cTMD(k1-nModes);
        CC2(k1,1)=-cTMD(k1-nModes);
    end
end  

CC = CC1 + CC2;
clear k1 k2


%% 导入气动力模型数据

for k1=1
ExpName=ExpNames(k1,:);
end

%ExpName:实验名称
%girderindex:梁号(上游输入1，下游输入2)
%TMDsindex:安装TMD的需要（TMD的序号，每个TMD对应的序号见TMD_logfile.mat）
%freshift:频率偏移，因为TMD频率分辨率很低，在一定范围中移动可以更好与试验结果相对照

% 导入试验数据
my_table=load('SZTD110_logfile.mat');
my_table=my_table.my_table;
fname = ExpName;
chanum = 8; %记录通道数
filename = strsplit(fname, '-');
casenumber = cell2mat(filename(3));
spacing = str2double(cell2mat(filename(4)));
type = cell2mat(filename(5));
Rspeed = cell2mat(filename(6));

casenumber = str2double(casenumber(5:end));
Rspeed = Rspeed(1:end - 1);
Rspeed = str2double(Rspeed);

if strcmp(type, 'fasan')
    type = 1;
else

    if strcmp(type, 'shuaijian')
        type = 2;
    else
        error("实验数据文件名可能有误，请确认信号状态为衰减或发散，程序终止")
    end

end

% 判断工况风攻角
if casenumber <= 17 || and(casenumber >= 32, casenumber <= 35) || and(casenumber >= 40, casenumber <= 16) || casenumber == 55
    AOA = 3;
else

    if and(casenumber >= 18, casenumber <= 25) || and(casenumber >= 36, casenumber <= 37) || and(casenumber >= 47, casenumber <= 50) || casenumber == 54
        AOA = 0;
    else

        if and(casenumber >= 26, casenumber <= 31) || and(casenumber >= 38, casenumber <= 39) || and(casenumber >= 51, casenumber <= 53) || and(casenumber >= 56, casenumber <= 57)
            AOA = -3;
        end

    end

end

if girderindex == 1
    girder = 'up';
    %导入试验数据
    Name = my_table.casename;
    isexist = find(Name == fname);
    a1 = my_table.up_parameter_a1(isexist);
    a2 = my_table.up_parameter_a2(isexist);
    a3 = my_table.up_parameter_a3(isexist);
    a4 = my_table.up_parameter_a4(isexist);
    a5 = my_table.up_parameter_a5(isexist);
    a = [a1 a2 a3 a4 a5];

    H4 = my_table.up_parameter_H4(isexist); % 气动刚度
    upperlimit = my_table.up_upperlimit(isexist); %除以特征长度D的无量纲振幅
    lowerlimit = my_table.up_lowerlimit(isexist); %除以特征长度D的无量纲振幅
    Fren_vibration_withwind = my_table.up_Fren_vibration_withwind(isexist);
    F0 = my_table.up_Fre_vibration(isexist); % Frequency without wind
    Zeta0 =  my_table.up_dltx_zeta0(isexist); % damping ratio without wind
    ReducedFrequency = my_table.up_Fre_vibration(isexist); % 折减频率
else
    girder = 'down';
    %导入试验数据
    Name = my_table.casename;
    isexist = find(Name == fname);
    a1 = my_table.down_parameter_a1(isexist);
    a2 = my_table.down_parameter_a2(isexist);
    a3 = my_table.down_parameter_a3(isexist);
    a4 = my_table.down_parameter_a4(isexist);
    a5 = my_table.down_parameter_a5(isexist);
    a = [a1 a2 a3 a4 a5];

    H4 = my_table.down_parameter_H4(isexist); % 气动刚度
    upperlimit = my_table.down_upperlimit(isexist); %除以特征长度D的无量纲振幅
    lowerlimit = my_table.down_lowerlimit(isexist); %除以特征长度D的无量纲振幅
    Fren_vibration_withwind = my_table.down_Fren_vibration_withwind(isexist);
    F0 = my_table.down_Fre_vibration(isexist); % Frequency without wind
    Zeta0 =  my_table.down_dltx_zeta0(isexist); % damping ratio without wind
    ReducedFrequency = my_table.down_Fre_vibration(isexist); % 折减频率
end

U=2*pi*Freq(mode_number)*D/ReducedFrequency; %风速
rho=1.225;
omega0=2*pi*F0; %无风振动频率
m=MM(mode_number,mode_number); %质量
b1=(rho*U*D*a1/m-2*Zeta0*omega0)*m;
b2=rho*U*a2;
b3=rho*U*a3/D;
b4=rho*U*a4/D^2;
b5=rho*U*a5/D^3;

%% 响应计算

h=0.1;% Time step
t=0:h:100;% Time
p=zeros(2,length(t)); %Initialize external load
gamma = 1/2; % Parameter in the Newmark algorithm
beta = 1/4; % Parameter in the Newmark algorithm
gfun = @(u,udot) bridge_damper(u,udot,rou,U,D,b1,b2,b3,b4,b5,MM,CC,KK,mode_integral_2,mode_integral_3,mode_integral_4,mode_integral_5,mode_integral_6,gamma,beta,h); % Handle to the nonlinear function
% TODO: 还没有考虑气动刚度

%% 所需要使用的函数
%% functions to be used
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [g,ks]=bridge_damper(u,udot,rou,U,D,a1,a2,a3,a4,a5,MM,CC,KK,mode_integral_2,mode_integral_3,mode_integral_4,mode_integral_5,mode_integral_6,gamma,beta,h,matrixsize,mode_number)
for k1 = 1:matrixsize
    g(k1)=0;
    if k1==mode_number
        for k2 = 1:matrixsize
            g(k1)=g(k1)+CC(k1,k2)*udot(k2)+KK(k1,k2)*u(k2);
        end
        g(k1)=g(k1)-b1*udot(k1)*
    
function result = P_eq(mode, temp_vec, Matrix)
    vec = temp_vec(:, mode);
    result = vec' * Matrix * vec;
end

function result = phiY(node, Mmapping, mode_vec, nModes)
    position_index = Mmapping.MatrixEqn(find(and(Mmapping.Node == node, Mmapping.DOF == 'UY')));

    if isempty(position_index)
        result = zeros(1, nModes);
    else
        result = mode_vec(position_index, 1:nModes);
    end

end

function result = Peq(Pmode, mode_vec, Mmapping, P_eachpoint, points, t)
    result = zeros(1, size(t, 2));

    for t1 = 1:points

        if sum(and(Mmapping.Node == t1, Mmapping.DOF == 'UY')) == 0
            result = result + 0 * P_eachpoint(t1, :);
        else
            position_index = Mmapping.MatrixEqn(find(and(Mmapping.Node == t1, Mmapping.DOF == 'UY')));
            result = result + mode_vec(position_index, Pmode) * P_eachpoint(t1, :);
        end

    end

end