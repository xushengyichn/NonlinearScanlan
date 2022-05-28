%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: Shengyi xushengyichn@outlook.com
%Date: 2022-05-22 11:22:16
%LastEditors: Shengyi xushengyichn@outlook.com
%LastEditTime: 2022-05-23 16:49:35
%FilePath: \NonlinearScanlan\NonlinearScanlan_Beam.m
%Description: 实现通过Scanlan气动力非线性模型计算简支梁的动力响应 Realization of the dynamic response of a simply supported beam calculated by the Scanlan aerodynamic nonlinear model
%
%Copyright (c) 2022 by Shengyi xushengyichn@outlook.com, All Rights Reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 检验瑞利阻尼矩阵中是否考虑TMD的刚度对结果的影响
% % Check whether the effect of TMD stiffness on the results is considered in the Rayleigh damping matrix
% Mian
% data = u(59, 1:round(end * 0.9, 0));
% save compare1 data
% Mian_FULLmatrix_unchange
% data2 = u(59, 1:round(end * 0.9, 0));
% close all
% load compare1
% clearvars -except data2 data
% plot(data(1:1000))
% hold on
% plot(data2(1:1000))
% title("Comparation of displacement with the Rayleigh damping matirx influence by the stiffness of TMD or not")
% legend("Regardless of TMD stiffness", "Regardless of TMD stiffness")
% % 结果表明影响非常大，但是对实际计算影响不大，实际的阻尼比需要实测获得
% % The results show that the influence is very large, but it has little influence on the actual calculation, and the actual damping ratio needs to be measured and obtained.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 设置TMD参数
%% Set TMD parameters
nTMD = 2;
mTMD = [100 100];
cTMD = [4.379642045362958 17.523840095542650];
kTMD = [19.181264445511033 307.0849728495260];
nodeTMD = [50 20]; %Node number(location of the TMD)
% mTMD = [2400];
% cTMD = [2 * mTMD(1) * 0.069704171453635 * 0.05];
% kTMD = [mTMD(1) * 0.069704171453635];
% nodeTMD = [50]; %Node number(location of the TMD)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 设置模态叠加法参数
%% Set the modal superposition method parameters
% Number of modes considered
nModes = 20; %考虑模态

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 设置运行参数Set run parameters

% 设置模态叠加法生成阻尼矩阵的方式
% 1 表示采用瑞利阻尼矩阵生成对应的模态叠加法阻尼矩阵
% 2 表示采用模态阻尼比直接生成的模态叠加法阻尼矩阵，需要指定模态阻尼比zeta
% Set the method for generating the damping matrix by the modal superposition method
% 1 means use the Rayleigh damping matrix to generate the corresponding modal superposition method damping matrix
% 2 means use the modal damping ratio to directly generate the modal superposition method damping matrix, you need to specify the modal Damping ratio zeta
DampingMatrixParameter = 1;

if DampingMatrixParameter == 2
    zeta = 0.3/100;
end

% 设置桥面坐标导入模式
% 1 采用程序自动写入一个含有桥面节点的数组，从 1 到 n
% 2 读入文件
% Set the bridge deck coordinate import mode
% 1 Use the program to automatically write an array containing the bridge deck nodes, from 1 to n
% 2 read in file
nodeondeckimport = 1;

if nodeondeckimport == 1
    %采用导入模式1时默认桥面节点为按顺序排列的的数列，只需要输入节点数
    % When importing mode 1 is adopted, the default bridge deck node is a sequence arranged in order, and only the number of nodes needs to be input.
    points = 101; % 主梁节点数 Number of girder nodes
end

% 设置外荷载施加方式
% 1 施加模态力
% 2 施加节点力（需要调整对应施加荷载部分的代码）
% Set the external load application method
% 1 Apply modal force
% 2 Apply nodal forces (need to adjust the code corresponding to the applied load section)
externalforcemethod = 2;

if externalforcemethod == 1
    % 施加模态力时需要定义模态力系数，模态力的频率为不施加TMD时的结构频率。
    % The modal force coefficient (MFC) needs to be defined when applying modal force, and the frequency of the modal force is the structural frequency when TMD is not applied.

    
    MFC=[1;0;0;0;0;0;0;0;0;0];
    if length(MFC)<nModes
        disp("荷载系数个数为: "+num2str(length(MFC))+"，小于模态数"+num2str(nModes)+"，其余模态力系数均补充为0")
        disp("原模态力系数: "+num2str(MFC'))
        temp_data=nModes-length(MFC);
        for t1 = 1:temp_data
            MFC = [MFC;0];
        end
        disp("现模态力系数: "+num2str(MFC'))
    end
end

% 是否输出振动视频
% Whether to output vibration video
% 1 export vides
% 2 do not export vides
exportvideo=2;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 将ANSYS中的稀疏矩阵处理为完全矩阵
%% Handling sparse matrices in ANSYS as full matrices
% 导入ANSYS MCK矩阵
% Import MCK matrix from ANSYS
hb_to_mm ('KMatrix.matrix', 'K.txt');
hb_to_mm ('MMatrix.matrix', 'M.txt');
hb_to_mm ('CMatrix.matrix', 'C.txt');
hb_to_mm ('KMatrixTMD.matrix', 'KTMD.txt');
%map the node and matrix from the KMatrix.mapping and MMatrix.mapping
Kdata = importdata('K.txt').data;
Kmatrix = zeros(Kdata(1, 1), Kdata(1, 2));

for i = 2:size(Kdata, 1)
    Kmatrix(Kdata(i, 1), Kdata(i, 2)) = Kdata(i, 3);
end

KTMDdata = importdata('KTMD.txt').data;
KTMDmatrix = zeros(KTMDdata(1, 1), KTMDdata(1, 2));

for i = 2:size(KTMDdata, 1)
    KTMDmatrix(KTMDdata(i, 1), KTMDdata(i, 2)) = KTMDdata(i, 3);
end

Mdata = importdata('M.txt').data;
Mmatrix = zeros(Mdata(1, 1), Mdata(1, 2));

for i = 2:size(Mdata, 1)
    Mmatrix(Mdata(i, 1), Mdata(i, 2)) = Mdata(i, 3);
end

Cdata = importdata('C.txt').data;
Cmatrix_DP = zeros(Cdata(1, 1), Cdata(1, 2));

for i = 2:size(Cdata, 1)
    Cmatrix_DP(Cdata(i, 1), Cdata(i, 2)) = Cdata(i, 3);
end

% 还原对角线以上元素，使之为对称阵, ANSYS只给出下三角矩阵
% Restore the elements above the diagonal to make it a symmetric matrix, ANSYS only gives the lower triangular matrix
K = diag(diag(Kmatrix) / 2) + Kmatrix - diag(diag(Kmatrix));
K = K + K';
M = diag(diag(Mmatrix) / 2) + Mmatrix - diag(diag(Mmatrix));
M = M + M';
C_exp = diag(diag(Cmatrix_DP) / 2) + Cmatrix_DP - diag(diag(Cmatrix_DP));
C_exp = C_exp + C_exp';
KTMD = diag(diag(KTMDmatrix) / 2) + KTMDmatrix - diag(diag(KTMDmatrix));
KTMD = KTMD + KTMD';
% 特征值分析，即计算频率Freq和振型Phi，calmodes数字代表求解的阶数，eigs中参数SM表示从较小的特征值开始求解
% Eigenvalue analysis, that is to calculate the frequency Freq and mode shape Phi, the calmodes number represents the order of the solution, and the parameter SM in eigs represents the solution from the smaller eigenvalue.
calmodes = 200; %考虑模态数 Consider the number of modes
[eig_vec, eig_val] = eigs(K, M, calmodes, 'SM');
[nfdof, nfdof] = size(eig_vec);

for j = 1:nfdof
    mnorm = sqrt(eig_vec(:, j)' * M * eig_vec(:, j));
    eig_vec(:, j) = eig_vec(:, j) / mnorm; %振型质量归一化 Mode shape mass normalization
end

[omeg, w_order] = sort(sqrt(diag(eig_val)));
mode_vec = eig_vec(:, w_order);
Freq = omeg / (2 * pi);

% 导入矩阵序号对应节点自由度关系。导入KCM中 *.mapping 文件的任意一个即可。他们是一样的。
% The imported matrix number corresponds to the node degree of freedom relationship. Import any one of the *.mapping files in KCM. They are the same.
KMmapping = importmappingmatrix('KMatrix.mapping');

% 仅保留需要的变量 M C K 矩阵 频率 圆频率 模态向量 对应关系
% Keep only the required variables M C K Matrix Frequency Circular frequency Mode vector Correspondence
clearvars -except M C K Freq omeg mode_vec C_exp KMmapping DampingMatrixParameter zeta nodeondeckimport points externalforcemethod nTMD nModes MFC mTMD cTMD kTMD nodeTMD exportvideo

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 设置TMD参数
%% Set TMD parameters
% mTMD = [100 100];
% cTMD = [2 * mTMD(1) * omeg(1) * 0.05 2 * mTMD(2) * omeg(2) * 0.05];
% kTMD = [mTMD(1) * omeg(1)^2 mTMD(2) * omeg(2)^2];
% nodeTMD = [50 20]; %Node number(location of the TMD)

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
KK = zeros(matrixsize, matrixsize);

for t1 = 1:matrixsize

    if t1 <= nModes
        KK(t1, t1) = P_eq(t1, mode_vec, K); %P=parameters
    end

    if and(t1 > nModes, t1 <= matrixsize)
        KK(t1, t1) = kTMD(t1 - nModes);
    end

end

clear t1

for t1 = 1:nModes

    for t2 = 1:nTMD

        for t3 = 1:nModes
            temp_data = kTMD(t2) * phiTMD(t2, t3) * phiTMD(t2, t1);
            KK(t1, t3) = KK(t1, t3) + temp_data;
            clear temp_data
        end

        indextmdt2 = nModes + t2;
        KK(t1, indextmdt2) = KK(t1, indextmdt2) - kTMD(t2) * phiTMD(t2, t1);
    end

end

clear t1 t2 t3

for t1 = 1:nTMD

    for t2 = 1:nModes
        temp_data = -kTMD(t1) * phiTMD(t1, t2);
        KK(nModes + t1, t2) = KK(nModes + t1, t2) + temp_data;
        clear temp_data
    end

end

clear t1

% 创建阻尼矩阵
% Create the Damping Matrix
CC = zeros(matrixsize, matrixsize);
CC = 0.013699749746335*KK;
switch DampingMatrixParameter
    case 1

        for t1 = 1:matrixsize

            if t1 <= nModes
%                 CC(t1, t1) = P_eq(t1, mode_vec, C_exp); %P=parameters
            end

            if and(t1 > nModes, t1 <= matrixsize)
                CC(t1, t1) =CC(t1, t1)+ cTMD(t1 - nModes);
            end

        end
        

        clear t1
    case 2

        for t1 = 1:matrixsize

            if t1 <= nModes
                CC(t1, t1) = 2 * MM(t1, t1) * omeg(t1) * zeta;
            end

            if and(t1 > nModes, t1 <= matrixsize)
                CC(t1, t1) = cTMD(t1 - nModes);
            end

        end

        clear t1
end

for t1 = 1:nModes

    for t2 = 1:nTMD

        for t3 = 1:nModes
            temp_data = cTMD(t2) * phiTMD(t2, t3) * phiTMD(t2, t1);
            CC(t1, t3) = CC(t1, t3) + temp_data;
            clear temp_data
        end

        indextmdt2 = nModes + t2;
        CC(t1, indextmdt2) = CC(t1, indextmdt2) - cTMD(t2) * phiTMD(t2, t1);
    end

end

clear t1 t2 t3

for t1 = 1:nTMD

    for t2 = 1:nModes
        temp_data = -cTMD(t1) * phiTMD(t1, t2);
        CC(nModes + t1, t2) = CC(nModes + t1, t2) + temp_data;
        clear temp_data
    end

end

clear t1

% 特征值分析（计算增加TMD后结构的模态和振型）
% Eigenvalue analysis (calculate the mode and mode shape of the structure after adding TMD)
clearvars -except KK MM CC matrixsize nModes nTMD mode_vec KMmapping nodeondeckimport points nodeTMD externalforcemethod MFC omeg exportvideo
calmodes = matrixsize; %考虑模态数 Consider the number of modes
[eig_vec, eig_val] = eigs(KK, MM, calmodes, 'SM');
[nfdof, nfdof] = size(eig_vec);

for k1 = 1:nfdof
    mnorm = sqrt(eig_vec(:, k1)' * MM * eig_vec(:, k1));
    eig_vec(:, k1) = eig_vec(:, k1) / mnorm; %振型质量归一化 Mode shape mass normalization
end

clear k1 ndof
[omeg2, w_order] = sort(sqrt(diag(eig_val)));
mode_vec2 = eig_vec(:, w_order);
Freq2 = omeg2 / (2 * pi);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 导入ansys增加TMD后的矩阵(这部分代码仅在前期作为验证代码使用，经过验证无误后不再使用)
%% Import the matrix after adding TMD in ansys (this part of the code is only used as the verification code in the early stage, and it will no longer be used after verification)

% hb_to_mm ( 'KMatrixTMD.matrix', 'KTMD.txt' );
% % hb_to_mm ( 'KMatrixTMDre.matrix', 'KTMD2.txt' );
% hb_to_mm ( 'MMatrixTMD.matrix', 'MTMD.txt' );
% hb_to_mm ( 'CMatrixTMD.matrix', 'CTMD.txt' );
% 
% % Kdata2 = importdata('KTMD2.txt').data;
% % Kmatrix2 = zeros(Kdata2(1,1),Kdata2(1,2));
% % for i = 2:size(Kdata2,1)
% %     Kmatrix2(Kdata2(i,1),Kdata2(i,2)) = Kdata2(i,3);
% % end
% % K2 =diag(diag(Kmatrix2)/2)+Kmatrix2-diag(diag(Kmatrix2));
% % K2 = K2+K2';
% 
% %%map the node and matrix from the KMatrix.mapping and MMatrix.mapping
% Kdata = importdata('KTMD.txt').data;
% Kmatrix = zeros(Kdata(1,1),Kdata(1,2));
% for i = 2:size(Kdata,1)
%     Kmatrix(Kdata(i,1),Kdata(i,2)) = Kdata(i,3);
% end
% 
% Mdata = importdata('MTMD.txt').data;
% Mmatrix = zeros(Mdata(1,1),Mdata(1,2));
% 
% for i = 2:size(Mdata,1)
%     Mmatrix(Mdata(i,1),Mdata(i,2)) = Mdata(i,3);
% end
% 
% Cdata = importdata('CTMD.txt').data;
% Cmatrix_DP = zeros(Cdata(1,1),Cdata(1,2));
% 
% for i = 2:size(Cdata,1)
%     Cmatrix_DP(Cdata(i,1),Cdata(i,2)) = Cdata(i,3);
% end
% % 还原对角线以上元素，使之为对称阵, ANSYS只给出下三角矩阵
% K =diag(diag(Kmatrix)/2)+Kmatrix-diag(diag(Kmatrix));
% K = K+K';
% M =diag(diag(Mmatrix)/2)+Mmatrix-diag(diag(Mmatrix));
% M = M+M';
% C =diag(diag(Cmatrix_DP)/2)+Cmatrix_DP-diag(diag(Cmatrix_DP));
% % C = C+C'+0.01*K2;
% 
% % 特征值分析，即计算频率Freq和振型Phi，50代表求50阶，SM表示从较小的特征值开始求解
% calmodes=matrixsize;%考虑模态
% [eig_vec,eig_val]=eigs(K,M,calmodes,'SM');
% [nfdof,nfdof]=size(eig_vec);
% for j=1:nfdof
%     mnorm=sqrt(eig_vec(:,j)'*M*eig_vec(:,j));
%     eig_vec(:,j)=eig_vec(:,j)/mnorm;%振型质量归一化
% end
% [omeg,w_order] =sort(sqrt(diag(eig_val)));
% mode_vec3=eig_vec(:,w_order);
% Freq3 = omeg/(2*pi);
% 
% matrixsize=nTMD+nModes;
% KMmapping3 = importmappingmatrix('KMatrixTMD.mapping');
% phiTMD3=zeros(nTMD,matrixsize);
% nodeTMD3=[10001 10002];
% % phiTMD row:TMD for each loaction column:the mode shape at the each
% % location of tmd
% for t1=1:nTMD
%     for t2=1:matrixsize
%         position_index=KMmapping3.MatrixEqn(find(and(KMmapping3.Node==nodeTMD3(t1),KMmapping3.DOF=='UY')));
%         phiTMD3(t1,t2)=mode_vec3(position_index,t2);
%     end
% end
% clear t1 t2

% V1=mode_vec;
% V2=mode_vec2;
% V3=mode_vec3;
% save vectest V1 V2 V3

% %调整ANSYS振型使之与matlab对应
% decidevalue2=zeros(size(phiTMD3,2),1);
% for t1=1:size(phiTMD3,2)
%     decidevalue=phiTMD3(1,t1)/mode_vec2(end-1,t1);
%     temp=abs(mode_vec2(end-1,t1))-abs(phiTMD3(1,t1));
%     if temp>10e-6
%         disp("模态"+num2str(t1)+"的完全矩阵和缩减矩阵tmd振型差异过大"+num2str(temp))
%     end
%     if decidevalue<0
% %         phiTMD3(:,t1)=phiTMD3(:,t1)*-1;
% %         mode_vec3(:,t1)=mode_vec3(:,t1)*-1;
%         mode_vec2(:,t1)=mode_vec2(:,t1)*-1;
%         decidevalue2(t1)=1;
%     end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate the response
%% 计算响应
nodeondeck = [];

switch nodeondeckimport
case 1
    nodeondeck = linspace(1, points, points);
case 2
    nodeondeck = readmatrix("nodeondeck.txt");
    points = length(nodeondeck); % 主梁节点数 Number of girder nodes
end

dt = 0.01;
% 计算时间（秒）Computation time (seconds)
T = 300;
NNT = T / dt;
t = 0:dt:T;
switch externalforcemethod
case 1 %作用模态力 
    PP = zeros(matrixsize, length(t)); % 生成模态力矩阵，全零 Generate modal force matrix, all zeros
    % TMD上不作用外力，因此荷载为零 No external force acts on the TMD, so the load is zero
    for t1 = 1:matrixsize

        if t1 <= nModes
            PP(t1, :) = MFC(t1)*sin(omeg(t1)*t);
        else
            PP(t1, :) = zeros(1, size(t, 2));
        end

    end    
    clear t1

case 2 %作用节点力
    fangdaxishu = 10;
    fangdaxishu2 = 100;
    zhenxing1 = 0.99950656036586;
    zhenxing2 = -0.92977648588903;
    PI = 3.14159265359;
    FRE = 0.069704171453635;
    FRE2 = 0.278900577315756;
    THETA = 2 * PI * FRE;
    THETA2 = 2 * PI * FRE2;
    
    P1 = fangdaxishu * zhenxing1 * sin(THETA * t);
    P2 = fangdaxishu2 * zhenxing2 * sin(THETA2 * t);
    P_eachpoint = zeros(points, length(t)); %生成外荷载时间序列矩阵，全零 Generate external loads time series matrix, all zeros
    P_eachpoint(50, :) = P1; % 修改对应节点处荷载Modify the load at the corresponding node
    P_eachpoint(20, :) = P2;
    PP = zeros(matrixsize, length(t)); % 生成模态力矩阵，全零 Generate modal force matrix, all zeros

    % 计算模态力 Calculate modal forces
    % TMD上不作用外力，因此荷载为零 No external force acts on the TMD, so the load is zero
    for t1 = 1:matrixsize

        if t1 <= nModes
            PP(t1, :) = Peq(t1, mode_vec, KMmapping, P_eachpoint, points, t);
        else
            PP(t1, :) = zeros(1, size(t, 2));
        end

    end
    clear t1
end




nbeta = 0.25;
ngam = 0.5;
u0 = zeros(matrixsize, 1);
udot0 = zeros(matrixsize, 1);
[u1 udot u2dot] = NewmarkInt(t, MM, CC, KK, PP, ngam, nbeta, u0, udot0);


Dis = zeros(points, length(t));

for t2 = 1:points
    pointnumber = nodeondeck(t2); %查看某个点的振动时程
    phiResult = phiY(pointnumber, KMmapping, mode_vec, nModes);
    Dis_temp = zeros(1, length(t));

    for t1 = 1:nModes
        Dis_temp = Dis_temp + phiResult(t1) .* u1(t1, :);
    end

    Dis(t2, :) = Dis_temp;
    clear Dis_temp t1
end

clear t2




figure
plot(t,Dis(50,:))
hold on
dataANSYS=readmatrix("T_DIS50.txt");
dataANSYS=dataANSYS(2:end,:);
plot(dataANSYS(1:round(end*0.9),1),dataANSYS(1:round(end*0.9),2),'r')
title("Displacement of the midpoint of the span")
legend("MATLAB","ANSYS")

figure
plot(t,u1(21,:))
hold on
dataANSYS=readmatrix("TMD1_DIS.txt");
dataANSYS=dataANSYS(2:end,:);
plot(dataANSYS(1:round(end*0.9),1),dataANSYS(1:round(end*0.9),2),'r')
title("Displacement of TMD1")
legend("MATLAB","ANSYS")

figure
plot(t,u1(22,:))
hold on
dataANSYS=readmatrix("TMD2_DIS.txt");
dataANSYS=dataANSYS(2:end,:);
plot(dataANSYS(1:round(end*0.9),1),dataANSYS(1:round(end*0.9),2),'r')
title("Displacement of TMD2")
legend("MATLAB","ANSYS")


%% 所需要使用的函数
%% functions to be used
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
