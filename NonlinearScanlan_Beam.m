%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: Shengyi xushengyichn@outlook.com
%Date: 2022-05-22 11:22:16
%LastEditors: Shengyi xushengyichn@outlook.com
%LastEditTime: 2022-05-23 13:57:50
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
% 设置运行参数Set run parameters

% 设置模态叠加法生成阻尼矩阵的方式
% 1 表示采用瑞利阻尼矩阵生成对应的模态叠加法阻尼矩阵
% 2 表示采用模态阻尼比直接生成的模态叠加法阻尼矩阵，需要指定模态阻尼比zeta
% Set the method for generating the damping matrix by the modal superposition method
% 1 means use the Rayleigh damping matrix to generate the corresponding modal superposition method damping matrix
% 2 means use the modal damping ratio to directly generate the modal superposition method damping matrix, you need to specify the modal Damping ratio zeta
DampingMatrixParameter = 2;

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 将ANSYS中的稀疏矩阵处理为完全矩阵
%% Handling sparse matrices in ANSYS as full matrices
% 导入ANSYS MCK矩阵
% Import MCK matrix from ANSYS
hb_to_mm ('KMatrix.matrix', 'K.txt');
hb_to_mm ('MMatrix.matrix', 'M.txt');
hb_to_mm ('CMatrix.matrix', 'C.txt');

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
clearvars -except M C K Freq omeg mode_vec C_exp KMmapping DampingMatrixParameter zeta nodeondeckimport points

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 设置TMD参数
%% Set TMD parameters
nTMD = 2;
mTMD = [100 100];
cTMD = [2 * mTMD(1) * omeg(1) * 0.05 2 * mTMD(2) * omeg(2) * 0.05];
kTMD = [mTMD(1) * omeg(1)^2 mTMD(2) * omeg(2)^2];
nodeTMD = [50 20]; %Node number(location of the TMD)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 考虑TMD振动响应的模态叠加法
%% Modal Superposition Method Considering TMD Vibration Response
% Number of modes considered
nModes = 10; %考虑模态
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

switch DampingMatrixParameter
    case 1

        for t1 = 1:matrixsize

            if t1 <= nModes
                CC(t1, t1) = P_eq(t1, mode_vec, C_exp); %P=parameters
            end

            if and(t1 > nModes, t1 <= matrixsize)
                CC(t1, t1) = cTMD(t1 - nModes);
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
clearvars -except KK MM CC matrixsize nModes nTMD mode_vec KMmapping nodeondeckimport points nodeTMD
calmodes = matrixsize; %考虑模态数 Consider the number of modes
[eig_vec, eig_val] = eigs(KK, MM, calmodes, 'SM');
[nfdof, nfdof] = size(eig_vec);

for k1 = 1:nfdof
    mnorm = sqrt(eig_vec(:, k1)' * MM * eig_vec(:, k1));
    eig_vec(:, k1) = eig_vec(:, k1) / mnorm; %振型质量归一化 Mode shape mass normalization
end

clear k1 ndof
[omeg, w_order] = sort(sqrt(diag(eig_val)));
mode_vec2 = eig_vec(:, w_order);
Freq2 = omeg / (2 * pi);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% 导入ansys增加TMD后的矩阵(这部分代码仅在前期作为验证代码使用，经过验证无误后不再使用)
% %% Import the matrix after adding TMD in ansys (this part of the code is only used as the verification code in the early stage, and it will no longer be used after verification)

% hb_to_mm ( 'KMatrixTMD.matrix', 'KTMD.txt' );
% hb_to_mm ( 'KMatrixTMDre.matrix', 'KTMD2.txt' );
% hb_to_mm ( 'MMatrixTMD.matrix', 'MTMD.txt' );
% hb_to_mm ( 'CMatrixTMD.matrix', 'CTMD.txt' );

% Kdata2 = importdata('KTMD2.txt').data;
% Kmatrix2 = zeros(Kdata2(1,1),Kdata2(1,2));
% for i = 2:size(Kdata2,1)
%     Kmatrix2(Kdata2(i,1),Kdata2(i,2)) = Kdata2(i,3);
% end
% K2 =diag(diag(Kmatrix2)/2)+Kmatrix2-diag(diag(Kmatrix2));
% K2 = K2+K2';

% %%map the node and matrix from the KMatrix.mapping and MMatrix.mapping
% Kdata = importdata('KTMD.txt').data;
% Kmatrix = zeros(Kdata(1,1),Kdata(1,2));
% for i = 2:size(Kdata,1)
%     Kmatrix(Kdata(i,1),Kdata(i,2)) = Kdata(i,3);
% end

% Mdata = importdata('MTMD.txt').data;
% Mmatrix = zeros(Mdata(1,1),Mdata(1,2));

% for i = 2:size(Mdata,1)
%     Mmatrix(Mdata(i,1),Mdata(i,2)) = Mdata(i,3);
% end

% Cdata = importdata('CTMD.txt').data;
% Cmatrix_DP = zeros(Cdata(1,1),Cdata(1,2));

% for i = 2:size(Cdata,1)
%     Cmatrix_DP(Cdata(i,1),Cdata(i,2)) = Cdata(i,3);
% end
% % 还原对角线以上元素，使之为对称阵, ANSYS只给出下三角矩阵
% K =diag(diag(Kmatrix)/2)+Kmatrix-diag(diag(Kmatrix));
% K = K+K';
% M =diag(diag(Mmatrix)/2)+Mmatrix-diag(diag(Mmatrix));
% M = M+M';
% C =diag(diag(Cmatrix_DP)/2)+Cmatrix_DP-diag(diag(Cmatrix_DP));
% C = C+C'+0.01*K2;

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

fangdaxishu = 10;
fangdaxishu2 = 100;
zhenxing1 = 0.99950656036586;
zhenxing2 = -0.92977648588903;
PI = 3.14159265359;
FRE = 0.069704171453635;
FRE2 = 0.278900577315756;
THETA = 2 * PI * FRE;
THETA2 = 2 * PI * FRE2;
dt = 0.1;
% 计算时间（秒）Computation time (seconds)
T = 500;
NNT = T / dt;
t = 0:dt:T;
P1 = fangdaxishu * zhenxing1 * sin(THETA * t)*10;
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

nbeta = 0.25;
ngam = 0.5;
u0 = zeros(matrixsize, 1);
udot0 = zeros(matrixsize, 1);
[u1 udot u2dot] = NewmarkInt(t, MM, CC, KK, PP, ngam, nbeta, u0, udot0);

tic
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
toc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% Create a new VideoWriter object (an empty video file). Use whatever format you want,
% but in my experience MP4 (with H.264 codec) is by far the best. Please stop using AVI.
hvid = VideoWriter('./movie.mp4', 'MPEG-4');

% Full quality, because why not?
set(hvid, 'Quality', 100);

% Set the frame rate
set(hvid, 'FrameRate', 30);

% Open the object for writing
open(hvid);

% Desired frame resolution (see fig2frame). The video will automatically adopt the resolution of the first frame (see HELP VIDEOWRITER).
% You could instead set the Width property of the video object, but I prefer this.
framepar.resolution = [1024, 768];

% Create a new figure
hfig = figure();
maxDis = max(max(abs(Dis)));

t_plot = downsample(t,10);
t_plot = t_plot(1:round(end/2,0));
for t1 = 1:length(t_plot)
    disp(['Processing frame ' num2str(t1) '...'])
    figure(hfig)
    t_temp = t_plot(t1);
    seq_temp = find(t==t_temp);
    Dis_temp = Dis(:, seq_temp) / maxDis;
    title("Displacement of the beam and TMD")
    plot(nodeondeck, Dis_temp)
    txt = ['Time: ' num2str(t_temp) ' s'];
    legend(txt)
    xlim([1 101])
    ylim([-1 1])
    hold on

    for t2 = 1:length(nodeTMD)
        Dis_tmd_temp = u1(nModes + t2, seq_temp);
        scatter(nodeTMD(t2), Dis_tmd_temp);
    end

    hold off
    % F = fig2frame(hfig,framepar); % <-- Use this
    F = getframe(hfig); % <-- Not this.
    clear Dis_temp
    writeVideo(hvid, F);
end

% Close the figure
close(hfig);
% Close the video object. This is important! The file may not play properly if you don't close it.
close(hvid);

% figure()
% plot(t,Dis,'b')
% xlabel('Time [s]'); ylabel("displacement of point:"+num2str(pointnumber));
% title("Displacement calculated by superposition method")

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
