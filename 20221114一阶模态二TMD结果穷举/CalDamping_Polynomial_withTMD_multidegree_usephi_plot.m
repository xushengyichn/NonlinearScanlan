%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%Date: 2022-09-26 19:35:04
%LastEditors: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%LastEditTime: 2022-10-19 12:16:19
%FilePath: \NonlinearScanlan\CalDamping_Polynomial_withTMD_multidegree.m
%Description: 计算多模态，施加某一阶模态多项式气动力模型后的响应，考虑TMD，通过干预大小振幅情况下的阻尼比，使得动力计算无需分段
%
%Copyright (c) 2022 by xushengyichn 54436848+xushengyichn@users.noreply.github.com, All Rights Reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clc; clear; close all;
function [result]=CalDamping_Polynomial_withTMD_multidegree_usephi_plot(nTMD,mTMD,zetaTMD,omegaTMD,modeTMD,mode_number,ifcalmode,MM_eq,KK_eq,calmodes,eig_val,eig_vec,u_temp)
% function modemaxdis_single=CalData_Polynomial_withTMD_multidegree(nTMD,mTMD,zetaTMD,omegaTMD,nodeTMD,mode_number,ifcalmode,calmodes,eig_val,eig_vec)
%% 参数设置

% nTMD; % TMD的数量   [1*1]
% mTMD; % TMD的质量 [1*n]
% zetaTMD; % TMD的阻尼比 [1*n]
% omegaTMD; % TMD的角频率 [1*n]
% nodeTMD; % TMD的节点 [1*n]
% mode_number; % 气动力施加的模态 [1*1]
% ifcalmode; % 是否计算模态 [1*1]     1:[脚本]：采用导入的KM矩阵计算模态； 2：[脚本]：直接导入计算好的特征值和特征向量； 3：[函数]：直接导入计算好的特征值和特征向量；
% MM_eq; % 原始模型模态化质量矩阵 [calmodes*calmodes]
% KK_eq; % 原始模型模态化刚度矩阵 [calmodes*calmodes]
% calmodes; % 计算模态的模态数 [1*1]
% eig_val; %原始模型的特征值
% eig_vec; %原始模型的特征向量
% t_length; %计算的时间长度


D = 20; %断面参考宽度

if ~exist("nTMD","var")
    nTMD = 0;
    disp("变量nTMD不存在，取默认值："+num2str(nTMD))
end
if ~exist("mTMD","var")
    mTMD = [2400] * 10;
    disp("变量mTMD不存在，取默认值："+num2str(mTMD))
end
if ~exist("zetaTMD","var")
    Ftmd = 0.833853594612216;
    omegatmd = 2 * pi * Ftmd;
    zetaTMD=0.05;
    cTMD = [2 * mTMD(1) * omegatmd * zetaTMD];
    disp("变量zetaTMD不存在，取默认值："+num2str(zetaTMD))
else
    cTMD = [2 * mTMD .* omegaTMD .* zetaTMD];
%     for k3 =1:length(mTMD)
%         cTMD(k3)=2 * mTMD(k3) .* omegaTMD(k3) .* zetaTMD(k3);
%     end
end
if ~exist("omegaTMD","var")
    Ftmd = 0.833853594612216;
    omegatmd = 2 * pi * Ftmd;
    kTMD = [mTMD(1) * omegatmd^2];
    disp("变量omegaTMD不存在，取默认值："+num2str(kTMD))
else
    kTMD = [mTMD.* omegaTMD.^2];
end
if ~exist("mode_number","var")
    mode_number = 1; %气动力施加的模态
    disp("变量mode_number不存在，取默认值："+num2str(mode_number))
end
% if ~exist("xTMD","var")
%     xTMD = [0]; %Node number(location of the TMD)
%     disp("变量xTMD不存在，取默认值："+num2str(xTMD))
% end


if ~exist("ifcalmode","var")
    ifcalmode = 1; %是否计算模态
    disp("变量ifcalmode不存在，取默认值："+num2str(ifcalmode))
end



if ifcalmode==1
    calmodes = 5; %考虑模态数 Consider the number of modes
    matrix=load('matrix.mat');
    K=matrix.K;
    M=matrix.M;
    disp("质量刚度矩阵不存在，自动导入中")
    [eig_vec, eig_val] = eigs(K, M, calmodes, 'SM');
else
    if ifcalmode==2
        calmodes = 5; %考虑模态数 Consider the number of modes
        modenew=load('modenew.mat');
        eig_val=modenew.eig_val;
        eig_vec=modenew.eig_vec;   
    else
%         disp("从函数读入计算完成的模态")
    end 
end

% mTMD = 0;
% cTMD = 0;
% kTMD = 0;
nModes = calmodes;

girderindex = 1;

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

Zeta0=0.3/100*ones(nModes,1);
% Zeta0(2) = 5/100;
% disp("模态阻尼比未定义，采用默认值，每阶模态阻尼比均为：" + num2str(Zeta0))



%% 读取梁桥模型

%% 将ANSYS中的稀疏矩阵处理为完全矩阵
%% Handling sparse matrices in ANSYS as full matrices
% 导入ANSYS MCK矩阵
% Import MCK matrix from ANSYS
% hb_to_mm ('KMatrix.matrix', 'K.txt');
% hb_to_mm ('MMatrix.matrix', 'M.txt');
% hb_to_mm ('CMatrix.matrix', 'C.txt');

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
exportvideo = 2;

% 特征值分析，即计算频率Freq和振型Phi，calmodes数字代表求解的阶数，eigs中参数SM表示从较小的特征值开始求解
% Eigenvalue analysis, that is to calculate the frequency Freq and mode shape Phi, the calmodes number represents the order of the solution, and the parameter SM in eigs represents the solution from the smaller eigenvalue.




[nfdof, nfdof] = size(eig_vec);

% for j = 1:nfdof
%     mnorm = sqrt(eig_vec(:, j)' * M * eig_vec(:, j));
%     eig_vec(:, j) = eig_vec(:, j) / mnorm; %振型质量归一化 Mode shape mass normalization
% end

% for j = 1:nfdof
%     m_modal(j) = sqrt(eig_vec(:, j)' * M * eig_vec(:, j));
% end


[omeg, w_order] = sort(sqrt(diag(eig_val)));
mode_vec = eig_vec(:, w_order);
Freq = omeg / (2 * pi);


%% 构建包含TMD的MCK矩阵
% 导入矩阵序号对应节点自由度关系。导入KCM中 *.mapping 文件的任意一个即可。他们是一样的。
% The imported matrix number corresponds to the node degree of freedom relationship. Import any one of the *.mapping files in KCM. They are the same.
KMmapping = importmappingmatrix('KMatrix.mapping');
UYNode=sort(KMmapping.Node(KMmapping.DOF == 'UY'));

% 计算桥面节点的振型向量
nodeondeck = importdata('nodeondeck.txt');
mode = zeros(length(nodeondeck), nModes);

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

for k1 = 1:length(nodegap) - 1
    nodedis(k1, 1) = nodegap(k1 + 1) - nodegap(k1);
end

clear k1
% 模态振型比节点间距多一个点，所以将两节点的模态取平均值计算振型积分
for k1 = 1:length(nodegap) - 1
    modecal(k1, :) = (mode(k1 + 1, :) + mode(k1, :)) / 2;
end

for k1 = 1:nModes
    integral_1(k1) = sum(abs(modecal(:, k1)) .* nodedis);
    integral_2(k1) = sum(modecal(:, k1).^2 .* nodedis);
    integral_3(k1) = sum(modecal(:, k1).^2 .* abs(modecal(:, k1)) .* nodedis);
    integral_4(k1) = sum(modecal(:, k1).^4 .* nodedis);
    integral_5(k1) = sum(modecal(:, k1).^2 .* abs(modecal(:, k1)).^3 .* nodedis);
    integral_6(k1) = sum(modecal(:, k1).^6 .* nodedis);
end

phi = mode(:, mode_number);
mode_integral_1 = integral_1(mode_number);
mode_integral_2 = integral_2(mode_number);
mode_integral_3 = integral_3(mode_number);
mode_integral_4 = integral_4(mode_number);
mode_integral_5 = integral_5(mode_number);
mode_integral_6 = integral_6(mode_number);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 考虑TMD振动响应的模态叠加法
%% Modal Superposition Method Considering TMD Vibration Response
% tic
matrixsize = nTMD + nModes;
phiTMD = modeTMD;
% phiTMD row:TMD for each loaction column:the mode shape at the each
% location of tmd
% 
% for t3 = 1:nTMD
% 
%     for t2 = 1:nModes
% 
%         [~,index]=sort(abs(nodegap-xTMD(t3)));%查找与xTMD最接近的点的排序
%         xResult=nodegap(index(1:2));%获取最接近的两个点的x坐标
%         mode2nodes=mode(index(1:2),1:nModes);%获取两个点坐标的y值
%         phi_result=interp1(xResult,mode2nodes,xTMD(t3),'linear','extrap');%插值以后任意点的振型
% %         disp(phi_result)
%         phiTMD(t3, t2) = phi_result(t2);
% 
%         
%         % % 找到TMD安装位置对应的振型向量大小并储存到phiTMD变量中
%         % % Find the size of the mode shape vector corresponding to the installation position of the TMD and store it in the phiTMD variable
%         % position_index = KMmapping.MatrixEqn(find(and(KMmapping.Node == nodeTMD(t1), KMmapping.DOF == 'UY')));
%         % if isempty(position_index)
%         %     phiTMD(t1, t2)=0;
%         % else
%         %     phiTMD(t1, t2) = mode_vec(position_index, t2);
%         % end
%     end
% 
% end
% clear t1
if nTMD ~= 0
%     disp("第一阶模态质量比" + num2str(phiTMD(1, 1)^2 * mTMD(1) * 100) + "%")
end

clear t1 t2

% 创建质量矩阵
% Create the Mass matrix
MM = zeros(matrixsize, matrixsize);
if ~exist("MM_eq","var")
    for t1 = 1:matrixsize

        if t1 <= nModes
            MM(t1, t1) = P_eq(t1, mode_vec, M); %行列数小于等于nModes为模态质量 The number of rows and columns is less than nModes is the modal quality
        end

        if and(t1 > nModes, t1 <= matrixsize)
            MM(t1, t1) = mTMD(t1 - nModes); %行列数大于nModes为TMD实际质量 The number of rows and columns is greater than nModes for the actual quality of TMD
        end

    end
else
    MM(1:nModes,1:nModes)=MM_eq;
    for k1 =nModes+1:matrixsize
        MM(k1,k1)=mTMD(k1-nModes);
    end
end
clear t1

for j = 1:nfdof
    m_modal(j) =MM(j,j);
end

% 创建刚度矩阵
% Create the Stiffness Matrix
KK = zeros(size(MM, 1), size(MM, 2));
KK1 = zeros(size(MM, 1), size(MM, 2));

if ~exist("KK_eq","var")
    for k1 = 1:matrixsize

        if k1 <= nModes
            KK1(k1, k1) = P_eq(k1, mode_vec, K); %P=parameters
        elseif k1 > nModes
            KK1(k1, k1) = kTMD(k1 - nModes);
        end

    end
else
    KK1(1:nModes,1:nModes)=KK_eq;
    for k1 =nModes+1:matrixsize
        KK1(k1,k1)=kTMD(k1-nModes);
    end
end

KK2 = zeros(matrixsize, matrixsize);

for k1 = 1:nModes % 第k1行

    for k2 = 1:nModes % 第k2列

        for k3 = 1:nTMD % 第k3个TMD
            KK2(k1, k2) = KK2(k1, k2) + kTMD(k3) * phiTMD(k3, k1) * phiTMD(k3, k2);
        end

    end

end

clear k1 k2

for k1 = 1:nModes

    for k2 = 1:nTMD
        KK2(k1, k2 + nModes) = -kTMD(k2) * phiTMD(k2, k1);
        KK2(k2 + nModes, k1) = -kTMD(k2) * phiTMD(k2, k1);
    end

end


KK = KK1 + KK2;
clear k1 k2

% 创建阻尼矩阵
% Create the Damping Matrix
CC = zeros(size(MM, 1), size(MM, 2));
CC1 = zeros(size(MM, 1), size(MM, 2));
CC2 = zeros(size(MM, 1), size(MM, 2));

for k1 = 1:matrixsize

    if k1 <= nModes
        CC1(k1, k1) = Zeta0(k1) * 4 * pi * MM(k1, k1) * Freq(k1);
    elseif k1 > nModes
        CC1(k1, k1) = cTMD(k1 - nModes);
    end

end

for k1 = 1:nModes % 第k1行

    for k2 = 1:nModes % 第k2列

        for k3 = 1:nTMD % 第k3个TMD
            CC2(k1, k2) = CC2(k1, k2) + cTMD(k3) * phiTMD(k3, k1) * phiTMD(k3, k2);
        end

    end

end

clear k1 k2

for k1 = 1:nModes

    for k2 = 1:nTMD
        CC2(k1, k2 + nModes) = -cTMD(k2) * phiTMD(k2, k1);
        CC2(k2 + nModes, k1) = -cTMD(k2) * phiTMD(k2, k1);
    end

end


CC = CC1 + CC2;
clear k1 k2
% toc

%% 导入气动力模型数据

for k1 = 1
    ExpName = ExpNames(k1, :);
end

%ExpName:实验名称
%girderindex:梁号(上游输入1，下游输入2)
%TMDsindex:安装TMD的需要（TMD的序号，每个TMD对应的序号见TMD_logfile.mat）
%freshift:频率偏移，因为TMD频率分辨率很低，在一定范围中移动可以更好与试验结果相对照

% 导入试验数据
my_table = load('SZTD110_logfile.mat');
my_table = my_table.my_table;
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
    % a = [a1 a2 a3 a4 a5];
    m_exp=my_table.up_m(isexist); %试验中单位长度模型的质量
    H4 = my_table.up_parameter_H4(isexist); % 气动刚度
    Fren_vibration_withwind = my_table.up_Fren_vibration_withwind(isexist);
    F0 = my_table.up_Fre_vibration(isexist); % Frequency without wind
    Zeta0 = my_table.up_dltx_zeta0(isexist); % damping ratio without wind
    ReducedFrequency = my_table.up_ReducedFre(isexist); % 折减频率
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
    % a = [a1 a2 a3 a4 a5];
    m_exp=my_table.down_m(isexist); %试验中单位长度模型的质量
    H4 = my_table.down_parameter_H4(isexist); % 气动刚度
    Fren_vibration_withwind = my_table.down_Fren_vibration_withwind(isexist);
    F0 = my_table.down_Fre_vibration(isexist); % Frequency without wind
    Zeta0 = my_table.down_dltx_zeta0(isexist); % damping ratio without wind
    ReducedFrequency = my_table.down_ReducedFre(isexist); % 折减频率
end

d1=m_modal(mode_number)*a1/m_exp;
d2=m_modal(mode_number)*a2/m_exp*phideckmax(mode_number);
d3=m_modal(mode_number)*a3/m_exp*phideckmax(mode_number)^2;
d4=m_modal(mode_number)*a4/m_exp*phideckmax(mode_number)^3;
d5=m_modal(mode_number)*a5/m_exp*phideckmax(mode_number)^4;





U = 2 * pi * Freq(mode_number) * D / ReducedFrequency; %风速
rho = 1.225;
omega0 = 2 * pi * Freq(mode_number); %无风振动频率
m = MM(mode_number, mode_number); %质量
b1 = rho * U * D * a1 / m_modal(mode_number);
b2 = rho * U * a2/ m_modal(mode_number);
b3 = rho * U * a3 / D/ m_modal(mode_number);
b4 = rho * U * a4 / D^2/ m_modal(mode_number);
b5 = rho * U * a5 / D^3/ m_modal(mode_number);

phi1max=max(mode_vec(:, 1));
m =  m_modal(mode_number);

u_temp=u_temp;
b_temp=(b1 * mode_integral_2 + b2 * abs(u_temp) * mode_integral_3 + b3 * u_temp^2 * mode_integral_4 + b4 * abs(u_temp)^3 * mode_integral_5 + b5 * u_temp^4 * mode_integral_6)/mode_integral_2;
t_temp_s=0:0.1:50;
for tt1= 1:length(t_temp_s)
b_temp_result(tt1)=(b1 * mode_integral_2 + b2 * abs(t_temp_s(tt1)) * mode_integral_3 + b3 * t_temp_s(tt1)^2 * mode_integral_4 + b4 * abs(t_temp_s(tt1))^3 * mode_integral_5 + b5 * t_temp_s(tt1)^4 * mode_integral_6)/mode_integral_2;
end
% figure
% plot(t_temp_s,b_temp_result)
a1_temp=b_temp/(rho * U * D)*m_modal(mode_number);
a1=a1_temp;




CC_aero=CC;
CC_aero(mode_number,mode_number)=CC_aero(mode_number,mode_number)-rho*U*D*a1*mode_integral_2;
KK(mode_number,mode_number)=KK(mode_number,mode_number)-rho*U*2*H4*mode_integral_2;%是否施加气动刚度，若不施加，注释本行
Mode_add_ADF = Complex_Eigenvalue_Analysis_complex(MM,CC_aero,KK);%增加线性气动力的模态
% disp(Mode)
Mode_sys = Complex_Eigenvalue_Analysis_complex(MM,CC,KK);%不考虑气动阻尼的阻尼比
% disp(Mode_sys)

result.Mode_add_ADF=Mode_add_ADF;
result.Mode_sys=Mode_sys;
% [eig_vec_mode, eig_val_mode] = eigs(KK, MM, matrixsize, 'SM');
% Fre_mode=eig_val_mode/2/pi;

function result = P_eq(mode, temp_vec, Matrix)
    vec = temp_vec(:, mode);
    result = vec' * Matrix * vec;
end



end
