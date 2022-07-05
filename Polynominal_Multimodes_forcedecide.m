%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: Shengyi xushengyichn@outlook.com
%Date: 2022-05-15 16:16:48
%LastEditors: Shengyi Xu xushengyichn@outlook.com
%LastEditTime: 2022-07-04 23:34:17
%FilePath: \NonlinearScanlan\Polynominal_Multimodes_forcedecide.m
%Description: 本函数目的为计算多自由度多项式模型响应，重点对比在荷载展开式中多阶模态对结果的影响
%
%Copyright (c) 2022 by Shengyi xushengyichn@outlook.com, All Rights Reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a function
% function out = NonlinearScanlan(Fre, Mass, Zeta0, rho, D, U, Y1k, epsilonk, Y2k, ClK, t, P, u0, udot0)

% Nonlinear Newmark's Direct Integration Method with Nonlinear Scanlan
% empirial noninear
% model++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% (n = number of time steps)
% (ndof = number degrees of freedom)

% INPUT
% Fre      = Frequency of the system         =>[1]
% Mass     = Mass of the system              =>[1]
% Zeta0    = Damping ratio of the system     =>[1]
% rho      = Air density                     =>[1]
% D        = Reference length                =>[1]
% U        = Wind speed at a certain reduced frequency       =>[1]
% Y1k      = Y1 at a certain reduced frequency       =>[1]
% epsilonk = Epsilon at a certain reduced frequency       =>[1]
% Y2k      = Y2 at a certain reduced frequency       =>[1]
% Clk      = Cl at a certain reduced frequency       =>[1]
% t        = Time vector         =>[1,n]
% P        = load vs. time       =>[1,n]
% u0       = Initial displacements =>[1]
% udot0    = Initial velocity =>[1]
% gam      = gamma (constant)
% beta     = beta  (constant)

%--------------------------------------------------------------------------
% beta = 0,     gamma = 1/2 => explicit central difference method
% beta = 1/4,   gamma = 1/2 => undamped trapezoidal rule (implicit)

%--------------------------------------------------------------------------
clc
clear
close all
i=0;
if i==1
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
%% 设置模态叠加法参数
%% Set the modal superposition method parameters
% Number of modes considered
nModes = 2; %考虑模态
mode = 1; %施加模态力
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 设置运行参数Set run parameters

% 设置模态叠加法生成阻尼矩阵的方式
% 1 表示采用瑞利阻尼矩阵生成对应的模态叠加法阻尼矩阵
% 2 表示采用模态阻尼比直接生成的模态叠加法阻尼矩阵，需要指定模态阻尼比zeta
% Set the method for generating the damping matrix by the modal superposition method
% 1 means use the Rayleigh damping matrix to generate the corresponding modal superposition method damping matrix
% 2 means use the modal damping ratio to directly generate the modal superposition method damping matrix, you need to specify the modal Damping ratio zeta
DampingMatrixParameter = 2;

if DampingMatrixParameter == 1
    beta = 0.013699749746335;
end

if DampingMatrixParameter == 2
    zeta = 0.3/100;
end

% 设置桥面坐标导入模式
% 1 采用程序自动写入一个含有桥面节点的数组，从 1 到 n
% 2 读入文件
% Set the bridge deck coordinate import mode
% 1 Use the program to automatically write an array containing the bridge deck nodes, from 1 to n
% 2 read in file
nodeondeckimport = 2;

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
externalforcemethod = 1;

if externalforcemethod == 1
    % 施加模态力时需要定义模态力系数，模态力的频率为不施加TMD时的结构频率。
    % The modal force coefficient (MFC) needs to be defined when applying modal force, and the frequency of the modal force is the structural frequency when TMD is not applied.

    MFC = [1; 0; 0; 0; 0; 0; 0; 0; 0; 0];

    if length(MFC) < nModes
        disp("荷载系数个数为: " + num2str(length(MFC)) + "，小于模态数" + num2str(nModes) + "，其余模态力系数均补充为0")
        disp("原模态力系数: " + num2str(MFC'))
        temp_data = nModes - length(MFC);

        for t1 = 1:temp_data
            MFC = [MFC; 0];
        end

        disp("现模态力系数: " + num2str(MFC'))
    end

end

% 是否输出振动视频
% Whether to output vibration video
% 1 export vides
% 2 do not export vides
exportvideo = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 将ANSYS中的稀疏矩阵处理为完全矩阵
%% Handling sparse matrices in ANSYS as full matrices
% 导入ANSYS MCK矩阵
% Import MCK matrix from ANSYS
hb_to_mm ('MMatrix.matrix', 'M.txt');
hb_to_mm ('KMatrix.matrix', 'K.txt');

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

% 还原对角线以上元素，使之为对称阵, ANSYS只给出下三角矩阵
% Restore the elements above the diagonal to make it a symmetric matrix, ANSYS only gives the lower triangular matrix
K = diag(diag(Kmatrix) / 2) + Kmatrix - diag(diag(Kmatrix));
K = K + K';
M = diag(diag(Mmatrix) / 2) + Mmatrix - diag(diag(Mmatrix));
M = M + M';

% 特征值分析，即计算频率Freq和振型Phi，calmodes数字代表求解的阶数，eigs中参数SM表示从较小的特征值开始求解
% Eigenvalue analysis, that is to calculate the frequency Freq and mode shape Phi, the calmodes number represents the order of the solution, and the parameter SM in eigs represents the solution from the smaller eigenvalue.
calmodes = 2; %考虑模态数 Consider the number of modes
[eig_vec, eig_val] = eigs(K, M, calmodes, 'SM');
[nfdof, nfdof] = size(eig_vec);

for j = 1:nfdof
    mnorm = sqrt(eig_vec(:, j)' * M * eig_vec(:, j));
    eig_vec(:, j) = eig_vec(:, j) / mnorm; %振型质量归一化 Mode shape mass normalization
end

[omeg, w_order] = sort(sqrt(diag(eig_val)));
mode_vec = eig_vec(:, w_order);
Freq = omeg / (2 * pi);

save temp
toc
end
load temp

% 导入矩阵序号对应节点自由度关系。导入KCM中 *.mapping 文件的任意一个即可。他们是一样的。
% The imported matrix number corresponds to the node degree of freedom relationship. Import any one of the *.mapping files in KCM. They are the same.
KMmapping = importmappingmatrix('KMatrix.mapping');

% 仅保留需要的变量 M C K 矩阵 频率 圆频率 模态向量 对应关系
% Keep only the required variables M C K Matrix Frequency Circular frequency Mode vector Correspondence
clearvars -except M C K Freq omeg mode_vec C_exp KMmapping DampingMatrixParameter zeta beta nodeondeckimport points externalforcemethod nTMD nModes MFC mTMD cTMD kTMD nodeTMD exportvideo mode

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 未施加TMD情况下响应的模态叠加法(仅考虑二阶，多项式模型仅采用最简单的Scanlan形式)
%% Modal Superposition Method for Responses Without TMD Applied
matrixsize = 2;

% 创建质量矩阵
% Create the Mass matrix
MM = zeros(matrixsize, matrixsize);

for t1 = 1:matrixsize
    MM(t1, t1) = P_eq(t1, mode_vec, M); %行列数小于等于nModes为模态质量 The number of rows and columns is less than nModes is the modal quality
end

clear t1

% 创建刚度矩阵
% Create the Stiffness Matrix
KK = zeros(matrixsize, matrixsize);

for t1 = 1:matrixsize
    KK(t1, t1) = P_eq(t1, mode_vec, K); %P=parameters
end

% 创建阻尼矩阵
% Create the Damping Matrix
CC = zeros(matrixsize, matrixsize);

for t1 = 1:matrixsize
    CC(t1, t1) = 2 * MM(t1, t1) * omeg(t1) * zeta;
end

b1 = 4e3;
b3 = -2e11;

gamma = 1/2; % Factor in the Newmark algorithm
beta = 1/4; % Factor in the Newmark algorithm



h = 0.01;
T = 100;
t = 0:h:T;
P = zeros(matrixsize, length(t));
pp = P;
% Y1 = Y1k;
% epsilon = epsilonk;
% Y2 = Y2k;

% MM = Mass;
% CC = Zeta0 * 4 * pi * MM * Fre;
% KK = 4 * pi^2 * MM * Fre^2;

%% Calculate the response
nodeondeck = [];

switch nodeondeckimport
    case 1
        nodeondeck = linspace(1, points, points);
    case 2
        nodeondeck = readmatrix("nodeondeck.txt");
        points = length(nodeondeck); % 主梁节点数 Number of girder nodes
        nodegap = readmatrix("nodegap.txt");
end

for t2 = 1:length(nodeondeck)
    pointnumber = nodeondeck(t2); %第n个节点
    phiResultall(:, t2) = phiY(pointnumber, KMmapping, mode_vec, nModes);
end

clear t2

k1=0;%是否进行不加tmd情况下的响应计算，1为是，0为否

if k1==1
% gfun = @(u, udot) Scanlan_nonlinear(u, udot, MM, CC, KK, gamma, beta, h, rho, U, D, Y1, epsilon, Y2, matrixsize, mode,phiResultall,nodegap);
gfun = @(u, udot) polysingle(u, udot, MM, CC, KK, gamma, beta, h, b1, b3, matrixsize, mode, phiResultall, nodegap);
u0 = zeros(matrixsize, 1);
u0_real = 0.01;
u0(mode) = u0_real / phiResultall(mode, 171);
u0(1) = 1e-1;
u0(2) = 1e-1;
% udot0 = zeros(matrixsize, 1);

    u = nonlinear_newmark_krenk(gfun, MM, pp, u0, udot0, gamma, beta, h);

    u_nounit = u(1, :);
    s = t * U / D;
    out = [s; u_nounit];
    figure
    plot(t, u(1, :))
    hold on
    plot(t, u(2, :))
    ylim([-1.5 1.5])
    title("Aerodynamic force only consider the first mode displacement")
    legend("Mode1","Mode2")


    gfun2 = @(u, udot) polymulti(u, udot, MM, CC, KK, gamma, beta, h, b1, b3, matrixsize, mode, phiResultall, nodegap);
    u0_real = 0.01;
    u0(mode) = u0_real / phiResultall(mode, 171);
    u0(1) = 1e-1;
    u0(2) = 1e-1;
    udot0 = zeros(matrixsize, 1);


    u = nonlinear_newmark_krenk(gfun2, MM, pp, u0, udot0, gamma, beta, h);

    u_nounit = u(1, :);
    s = t * U / D;
    out = [s; u_nounit];
    figure
    plot(t, u(1, :))
    hold on
    plot(t, u(2, :))
    ylim([-1.5 1.5])
    title("Aerodynamic force consider the first and second mode displacement")
    legend("Mode1","Mode2")
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 设置TMD参数
%% Set TMD parameters
nTMD = 1;
mTMD = [40000];
cTMD = [2 * mTMD(1) * 5.239256655795033 * 0.05];
kTMD = [mTMD(1) * 5.239256655795033^2];
nodeTMD = [1168]; %Node number(location of the TMD)

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

mTMD*phiTMD(1)^2

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
        CC = beta * KK;

        for t1 = 1:matrixsize

            if t1 <= nModes
                % CC(t1, t1) = P_eq(t1, mode_vec, C_exp); %P=parameters
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

% 考虑施加TMD的情况
P = zeros(matrixsize, length(t));
pp = P;
gfun = @(u, udot) polymulti_tmd(u, udot, MM, CC, KK, gamma, beta, h, b1, b3, matrixsize, mode, phiResultall, nodegap);
u0 = zeros(matrixsize, 1);
u0_real = 0.01;
u0(mode) = u0_real / phiResultall(mode, 171);
u0(1) = 1e-1;
u0(2) = 1e-1;
u0(3) = 0;
udot0 = zeros(matrixsize, 1);

u = nonlinear_newmark_krenk(gfun, MM, pp, u0, udot0, gamma, beta, h);




figure
plot(t, u(1, :))
hold on
plot(t, u(2, :))
plot(t, u(3,:))
ylim([-1.5 1.5])
title("Aerodynamic force only consider the first and second modes of displacement")
legend("Mode1","Mode2")
% [psd_avg, f, psd_plot] = fft_transfer(100,u(1, :))
% plot(f,psd_plot)

% % figure
% % [psd_avg, f, psd_plot] = fft_transfer(1/h,u_nounit');
% % plot(f,psd_plot)

% Dis = zeros(length(nodeondeck), length(t));

% for t2 = 1:length(nodeondeck)
%     pointnumber = nodeondeck(t2); %查看某个点的振动时程
%     phiResult = phiY(pointnumber, KMmapping, mode_vec, nModes);
%     Dis_temp = zeros(1, length(t));

%     for t1 = 1:nModes
%         Dis_temp = Dis_temp + phiResult(t1) .* u(t1, :);
%     end

%     Dis(t2, :) = Dis_temp;

%     clear Dis_temp t1
% end

% clear t2

% figure
% plot(t,Dis(171,:))

% [psd_avg, f, psd_plot]=fft_transfer(1/h,Dis(100,end/2:end)');
% figure
% plot(f, psd_plot)

% if nTMD>0
%     figure()
% plot(t,u(end,:))

% phitmd=phiY(nodeTMD, KMmapping, mode_vec, 1);
% disp("TMD的质量比为"+num2str(mTMD*phitmd^2))
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure
% % plot(t(1:end/10),Dis(171,1:end/10))
% plot(t,Dis(171,:))
% title("Displacement of node 171")
%     figure()
% plot(t,u(end,:))
% % plot(t,u(end,1:end/10))
% title("Displacement of TMD, mass ratio: "+num2str(mTMD*phitmd^2*100)+"%")

% if exportvideo ==1
% % Create a new VideoWriter object (an empty video file). Use whatever format you want,
% % but in my experience MP4 (with H.264 codec) is by far the best. Please stop using AVI.
% hvid = VideoWriter('./movie.mp4', 'MPEG-4');

% % Full quality, because why not?
% set(hvid, 'Quality', 100);

% % Set the frame rate
% set(hvid, 'FrameRate', 30);

% % Open the object for writing
% open(hvid);

% % Desired frame resolution (see fig2frame). The video will automatically adopt the resolution of the first frame (see HELP VIDEOWRITER).
% % You could instead set the Width property of the video object, but I prefer this.
% framepar.resolution = [1024, 768];

% % Create a new figure
% hfig = figure();
% maxDis = max(max(abs(Dis)));

% t_plot = downsample(t,10);
% t_plot = t_plot(1:round(end/2,0));
% for t1 = 1:length(t_plot)
%     disp(['Processing frame ' num2str(t1) '...'])
%     figure(hfig)
%     t_temp = t_plot(t1);
%     seq_temp = find(t==t_temp);
%     Dis_temp = Dis(:, seq_temp) / maxDis;
%     title("Displacement of the beam and TMD")
%     plot(nodeondeck, Dis_temp)
%     txt = ['Time: ' num2str(t_temp) ' s'];
%     legend(txt)
%     xlim([nodeondeck(1) nodeondeck(end)])
%     ylim([-1 1])
%     hold on

%     if nTMD > 0
%         for t2 = 1:length(nodeTMD)
%             Dis_tmd_temp = u(nModes + t2, seq_temp);
%             scatter(nodeTMD(t2), Dis_tmd_temp);
%         end
%         clear t2
%     end

%     hold off
%     % F = fig2frame(hfig,framepar); % <-- Use this
%     F = getframe(hfig); % <-- Not this.
%     clear Dis_temp
%     writeVideo(hvid, F);
% end
% clear t1

% % Close the figure
% close(hfig);
% % Close the video object. This is important! The file may not play properly if you don't close it.
% close(hvid);
% end

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

%% Nonlinear Newmark algorithm
function [u, udot, u2dot] = nonlinear_newmark_krenk(gfun, MM, pp, u0, udot0, gamma, beta, h)
    % Initialize variables
    u = zeros(size(MM, 1), size(pp, 2));
    udot = zeros(size(MM, 1), size(pp, 2));
    u2dot = zeros(size(MM, 1), size(pp, 2));
    % Initial conditions
    u(:, 1) = u0; % Initial displacements
    udot(:, 1) = udot0; % Initial velocities
    g = gfun(u0, udot0); % Evaluate the nonlinear function
    u2dot(:, 1) = MM \ (pp(:, 1) - g); % Calculate the initial accelerations

    for ii = 1:size(pp, 2) - 1
        % Prediction step
        u2dot(:, ii + 1) = u2dot(:, ii); % Predicted accelerations
        udot(:, ii + 1) = udot(:, ii) + h * u2dot(:, ii); % Predicted velocities
        u(:, ii + 1) = u(:, ii) + h * udot(:, ii) + 1/2 * h^2 * u2dot(:, ii); % Predicted displacements
        Nit = 0; % Number of iterations
        konv = 0; % Has the iterations converged 0=No, 1=Yes
        % Reusidal calculation, system matrices and increment correction
        while Nit < 1000 && konv == 0
            Nit = Nit + 1;
            [g, Ks] = gfun(u(:, ii + 1), udot(:, ii + 1)); % Calculate function value and the tangent
            rr = pp(:, ii + 1) - MM * u2dot(:, ii + 1) - g; % Calculate residual
            du = Ks \ rr; % Increment correction
            u(:, ii + 1) = u(:, ii + 1) + du; % Add incremental correction to the displacement
            udot(:, ii + 1) = udot(:, ii + 1) + gamma * h / (beta * h^2) * du; % Incremental correction for velocities
            u2dot(:, ii + 1) = u2dot(:, ii + 1) + 1 / (beta * h^2) * du; % Incremental correction accelerations

            if sqrt(rr' * rr) / length(rr) < 1.0e-4 % Convergence criteria
                konv = 1; % konv = 1 if the convergence criteria is fulfilled.
            end

        end

    end

end

%% Function file for the Scanlan nonlinear model
% function [g, ks] = Scanlan_nonlinear(u, udot, MM, CC, KK, gamma, beta, h, rho, U, D, Y1, epsilon, Y2, matrixsize, mode, phiResultall, nodegap)
%     Y1_diagMatrix = zeros(matrixsize, matrixsize); %Aerodynamice Force decision matrix 决定哪些自由度是否施加气动力
%     Y2_diagMatrix = zeros(matrixsize, matrixsize); %Aerodynamice Force decision matrix 决定哪些自由度是否施加气动力
%     epsilon_diagMatrix = zeros(matrixsize, matrixsize); %Aerodynamice Force decision matrix 决定哪些自由度是否施加气动力
%     diagMatrix = zeros(matrixsize, matrixsize); %Aerodynamice Force decision matrix 决定哪些自由度是否施加气动力
%     Y1_diagMatrix(mode, mode) = Y1(mode);
%     Y2_diagMatrix(mode, mode) = Y2(mode);
%     epsilon_diagMatrix(mode, mode) = epsilon(mode);
%     diagMatrix(mode, mode) = 1;

%     if ~exist('phiResultall', 'var')
%         phi2 = 1;
%         phi4 = 1;
%         disp("未输入振型，为单自由度问题。")
%     else
%         phi = phiResultall(mode, :);

%         for k1 = 1:length(nodegap) - 1
%             dx(k1) = nodegap(k1 + 1) - nodegap(k1);
%         end

%         clear k1
%         phi2 = 0;
%         phi4 = 0;

%         for k1 = 1:length(dx)
%             phi2 = phi2 + phi(k1)^2 * dx(k1);
%             phi4 = phi4 + phi(k1)^4 * dx(k1);
%         end

%         clear k1

%     end

%     g = (CC - rho .* U .* D .* Y1_diagMatrix .* (phi2 .* diagMatrix - epsilon_diagMatrix .* u.^2 .* phi4 ./ D.^2)) * udot + (KK) * u; % Function value
%     kc = CC - rho .* U .* D .* Y1_diagMatrix .* (phi2 .* diagMatrix - epsilon_diagMatrix .* u.^2 .* phi4 ./ D.^2);
%     u_udot = zeros(matrixsize, matrixsize);

%     for k1 = 1:size(u_udot, 1)
%         u_udot(k1, k1) = u(k1) * udot(k1);
%     end

%     kspring = KK + 2 .* rho .* U .* D .* Y1_diagMatrix .* epsilon_diagMatrix .* phi4 ./ D^2 * u_udot;
%     ks = kspring + gamma .* h ./ (beta .* h.^2) .* kc + 1 ./ (beta .* h.^2) .* MM; % Linearization

% end

%% Function file for the Scanlan nonlinear model
function [g, ks] = polysingle(u, udot, MM, CC, KK, gamma, beta, h, b1, b3, matrixsize, mode, phiResultall, nodegap)
    % Y1_diagMatrix = zeros(matrixsize, matrixsize);%Aerodynamice Force decision matrix 决定哪些自由度是否施加气动力
    % Y2_diagMatrix = zeros(matrixsize, matrixsize);%Aerodynamice Force decision matrix 决定哪些自由度是否施加气动力
    % epsilon_diagMatrix = zeros(matrixsize, matrixsize);%Aerodynamice Force decision matrix 决定哪些自由度是否施加气动力
    % diagMatrix = zeros(matrixsize, matrixsize);%Aerodynamice Force decision matrix 决定哪些自由度是否施加气动力

    phi1 = phiResultall(1, :);

    for k1 = 1:length(nodegap) - 1
        dx(k1) = nodegap(k1 + 1) - nodegap(k1);
    end

    clear k1
    phi1_2 = 0;
    phi1_4 = 0;

    for k1 = 1:length(dx)
        phi1_2 = phi1_2 + phi1(k1)^2 * dx(k1);
        phi1_4 = phi1_4 + phi1(k1)^4 * dx(k1);
    end

    clear k1

    q1 = u(1);
    q2 = u(2);
    qdot1 = udot(1);
    qdot2 = udot(2);

    g1 = CC(1, 1) * qdot1 - qdot1*b3*phi1_4*q1^2 -qdot1*b1*phi1_2 + KK(1, 1) * q1;
    g2 = CC(2, 2) * qdot2 + KK(2, 2) * q2;
    g = [g1; g2];


    kc1=-b3*q1^2*phi1_4-b1*phi1_2;
    kc2=0;
    kc3=0;
    kc4=0;
    kc = CC + [kc1 kc2 ;kc3 kc4 ];

%     kc = CC - [b3 * phi1_4 * q1^2+b1*phi1_2 0; 0 0];
    kspring = KK - [2 * b3 * phi1_4 * q1 * qdot1 0; 0 0];
    ks = kspring + gamma .* h ./ (beta .* h.^2) .* kc + 1 ./ (beta .* h.^2) .* MM; % Linearization

end

function [g, ks] = polymulti(u, udot, MM, CC, KK, gamma, beta, h, b1, b3, matrixsize, mode, phiResultall, nodegap)
    % Y1_diagMatrix = zeros(matrixsize, matrixsize);%Aerodynamice Force decision matrix 决定哪些自由度是否施加气动力
    % Y2_diagMatrix = zeros(matrixsize, matrixsize);%Aerodynamice Force decision matrix 决定哪些自由度是否施加气动力
    % epsilon_diagMatrix = zeros(matrixsize, matrixsize);%Aerodynamice Force decision matrix 决定哪些自由度是否施加气动力
    % diagMatrix = zeros(matrixsize, matrixsize);%Aerodynamice Force decision matrix 决定哪些自由度是否施加气动力

    phi1 = phiResultall(1, :);
    phi2= phiResultall(2, :);
    phi1= zeros(1,length(phi2));


    for k1 = 1:length(nodegap) - 1
        dx(k1) = nodegap(k1 + 1) - nodegap(k1);
    end

    clear k1
    phi1_1=0;
    phi2_1=0;
    phi1_2 = 0;
    phi2_2= 0;
    phi1_4=0;
    phi2_4=0;
    phi1_3_phi2=0;
    phi1_2_phi2_2=0;
    phi2_3_phi1=0;
    phi1_phi2=0;
    % phi1_3=0;
    % phi2_3=0;

    for k1 = 1:length(dx)
        phi1_1= phi1_1+phi1(k1)^1 * dx(k1);
        phi2_1= phi2_1+phi2(k1)^1 * dx(k1);
        phi1_2 = phi1_2 + phi1(k1)^2 * dx(k1);
        phi1_4 = phi1_4 + phi1(k1)^4 * dx(k1);
        phi2_2 = phi2_2 + phi2(k1)^2 * dx(k1);
        phi2_4 = phi2_4 + phi2(k1)^4 * dx(k1);
        % phi1_3 = phi1_3 + phi1(k1)^3 * dx(k1);
        % phi2_3 = phi2_3 + phi2(k1)^3 * dx(k1);
        phi1_3_phi2=phi1_3_phi2+phi1(k1)^3*phi2(k1)*dx(k1);
        phi1_2_phi2_2=phi1_2_phi2_2+phi1(k1)^2*phi2(k1)^2*dx(k1);
        phi2_3_phi1=phi2_3_phi1+phi2(k1)^3*phi1(k1)*dx(k1);
        phi1_phi2=phi1_phi2+phi1(k1)*phi2(k1)*dx(k1);
    end

    clear k1

    q1 = u(1);
    q2 = u(2);
    qdot1 = udot(1);
    qdot2 = udot(2);

    g1 = -qdot1*b3*phi1_4*q1^2 -2*qdot1*b3*phi1_3_phi2*q1*q2-qdot1*b3*phi1_2_phi2_2*q2^2-qdot2*b3*phi1_3_phi2*q1^2-2*qdot2*b3*phi1_2_phi2_2*q1*q2-qdot2*b3*phi2_3_phi1*q2^2-qdot1*b1*phi1_2-qdot2*b1*phi1_phi2+CC(1,1)*qdot1+KK(1,1)*q1; 
    g2 =-qdot2*b3*phi1_2_phi2_2*q1^2-2*qdot2*b3*phi2_3_phi1*q1*q2-qdot2*b3*phi2_4*q2^2-qdot2*b1*phi2_2+CC(2,2)*qdot2-qdot1*b3*phi1_3_phi2*q1^2-2*qdot1*b3*phi1_2_phi2_2*q1*q2-qdot1*b3*phi2_3_phi1*q2^2-qdot1*b1*phi1_phi2+KK(2,2)*q2 ;
    g = [g1; g2];


    kc1=-b3*q1^2*phi1_4-2*b3*q1*q2*phi1_3_phi2-b3*phi1_2_phi2_2*q2^2-b1*phi1_2;
    kc2=-b3*phi1_3_phi2*q1^2-2*b3*phi1_2_phi2_2*q1*q2-b3*phi2_3_phi1*q2^2-b1*phi1_phi2;
    kc3=-b3*phi1_3_phi2*q1^2-2*b3*phi1_2_phi2_2*q1*q2-b3*phi2_3_phi1*q2^2-b1*phi1_phi2;
    kc4=-b3*phi1_2_phi2_2*q1^2-2*b3*q1*q2*phi2_3_phi1-b3*q2^2*phi2_4-b1*phi2_2;
    kc = CC + [kc1 kc2 ;kc3 kc4 ];

    ks1=-2*b3*phi1_4*q1*qdot1-2*b3*phi1_3_phi2*qdot1*q2-2*b3*phi1_3_phi2*qdot2*q1-2*b3*phi1_2_phi2_2*qdot2*q2;
    ks2=-2*b3*phi1_3_phi2*qdot1*q1-2*b3*phi1_2_phi2_2*qdot1*q2-2*b3*phi1_2_phi2_2*qdot2*q1-2*b3*phi2_3_phi1*qdot2*q2;
    ks3=-2*b3*phi1_3_phi2*qdot1*q1-2*b3*phi1_2_phi2_2*qdot1*q2-2*b3*phi1_2_phi2_2*qdot2*q1-2*b3*phi2_3_phi1*qdot2*q2;
    ks4=-2*b3*phi1_2_phi2_2*qdot1*q1-2*b3*phi2_3_phi1*qdot1*q2-2*b3*phi2_3_phi1*qdot2*q1-2*b3*phi2_4*qdot2*q2;

    kspring = KK + [ks1 ks2;ks3 ks4];
    ks = kspring + gamma .* h ./ (beta .* h.^2) .* kc + 1 ./ (beta .* h.^2) .* MM; % Linearization

end


function [g, ks] = polymulti_tmd(u, udot, MM, CC, KK, gamma, beta, h, b1, b3, matrixsize, mode, phiResultall, nodegap)
    % Y1_diagMatrix = zeros(matrixsize, matrixsize);%Aerodynamice Force decision matrix 决定哪些自由度是否施加气动力
    % Y2_diagMatrix = zeros(matrixsize, matrixsize);%Aerodynamice Force decision matrix 决定哪些自由度是否施加气动力
    % epsilon_diagMatrix = zeros(matrixsize, matrixsize);%Aerodynamice Force decision matrix 决定哪些自由度是否施加气动力
    % diagMatrix = zeros(matrixsize, matrixsize);%Aerodynamice Force decision matrix 决定哪些自由度是否施加气动力

    phi1 = phiResultall(1, :);
    phi2 = phiResultall(2, :);
    phi2= zeros(1,length(phi2));


    for k1 = 1:length(nodegap) - 1
        dx(k1) = nodegap(k1 + 1) - nodegap(k1);
    end

    clear k1
    phi1_1=0;
    phi2_1=0;
    phi1_2 = 0;
    phi2_2= 0;
    phi1_4=0;
    phi2_4=0;
    phi1_3_phi2=0;
    phi1_2_phi2_2=0;
    phi2_3_phi1=0;
    phi1_phi2=0;
    % phi1_3=0;
    % phi2_3=0;

    for k1 = 1:length(dx)
        phi1_1= phi1_1+phi1(k1)^1 * dx(k1);
        phi2_1= phi2_1+phi2(k1)^1 * dx(k1);
        phi1_2 = phi1_2 + phi1(k1)^2 * dx(k1);
        phi1_4 = phi1_4 + phi1(k1)^4 * dx(k1);
        phi2_2 = phi2_2 + phi2(k1)^2 * dx(k1);
        phi2_4 = phi2_4 + phi2(k1)^4 * dx(k1);
        % phi1_3 = phi1_3 + phi1(k1)^3 * dx(k1);
        % phi2_3 = phi2_3 + phi2(k1)^3 * dx(k1);
        phi1_3_phi2=phi1_3_phi2+phi1(k1)^3*phi2(k1)*dx(k1);
        phi1_2_phi2_2=phi1_2_phi2_2+phi1(k1)^2*phi2(k1)^2*dx(k1);
        phi2_3_phi1=phi2_3_phi1+phi2(k1)^3*phi1(k1)*dx(k1);
        phi1_phi2=phi1_phi2+phi1(k1)*phi2(k1)*dx(k1);
    end

    clear k1

    q1 = u(1);
    q2 = u(2);
    utmd=u(3);
    qdot1 = udot(1);
    qdot2 = udot(2);
    udottmd=udot(3);

    g1 = -qdot1*b3*phi1_4*q1^2 -2*qdot1*b3*phi1_3_phi2*q1*q2-qdot1*b3*phi1_2_phi2_2*q2^2-qdot2*b3*phi1_3_phi2*q1^2-2*qdot2*b3*phi1_2_phi2_2*q1*q2-qdot2*b3*phi2_3_phi1*q2^2-qdot1*b1*phi1_2-qdot2*b1*phi1_phi2+CC(1,1)*qdot1+CC(1,2)*qdot2+CC(1,3)*udottmd+KK(1,1)*q1+KK(1,2)*q2+KK(1,3)*utmd; 
    g2 =-qdot2*b3*phi1_2_phi2_2*q1^2-2*qdot2*b3*phi2_3_phi1*q1*q2-qdot2*b3*phi2_4*q2^2-qdot2*b1*phi2_2+CC(2,2)*qdot2-qdot1*b3*phi1_3_phi2*q1^2-2*qdot1*b3*phi1_2_phi2_2*q1*q2-qdot1*b3*phi2_3_phi1*q2^2-qdot1*b1*phi1_phi2+KK(2,2)*q2 +CC(2,1)*qdot1+CC(2,3)*udottmd+KK(2,1)*q1+KK(2,3)*utmd;
    g3=CC(3,1)*qdot1+CC(3,2)*qdot2+CC(3,3)*udottmd+KK(3,1)*q1+KK(3,2)*q2+KK(3,3)*utmd;
    g = [g1; g2;g3];


    kc1=-b3*q1^2*phi1_4-2*b3*q1*q2*phi1_3_phi2-b3*phi1_2_phi2_2*q2^2-b1*phi1_2;
    kc2=-b3*phi1_3_phi2*q1^2-2*b3*phi1_2_phi2_2*q1*q2-b3*phi2_3_phi1*q2^2-b1*phi1_phi2;
    kc3=-b3*phi1_3_phi2*q1^2-2*b3*phi1_2_phi2_2*q1*q2-b3*phi2_3_phi1*q2^2-b1*phi1_phi2;
    kc4=-b3*phi1_2_phi2_2*q1^2-2*b3*q1*q2*phi2_3_phi1-b3*q2^2*phi2_4-b1*phi2_2;
    kc = CC + [kc1 kc2 0;kc3 kc4 0; 0 0 0];

    ks1=-2*b3*phi1_4*q1*qdot1-2*b3*phi1_3_phi2*qdot1*q2-2*b3*phi1_3_phi2*qdot2*q1-2*b3*phi1_2_phi2_2*qdot2*q2;
    ks2=-2*b3*phi1_3_phi2*qdot1*q1-2*b3*phi1_2_phi2_2*qdot1*q2-2*b3*phi1_2_phi2_2*qdot2*q1-2*b3*phi2_3_phi1*qdot2*q2;
    ks3=-2*b3*phi1_3_phi2*qdot1*q1-2*b3*phi1_2_phi2_2*qdot1*q2-2*b3*phi1_2_phi2_2*qdot2*q1-2*b3*phi2_3_phi1*qdot2*q2;
    ks4=-2*b3*phi1_2_phi2_2*qdot1*q1-2*b3*phi2_3_phi1*qdot1*q2-2*b3*phi2_3_phi1*qdot2*q1-2*b3*phi2_4*qdot2*q2;

    kspring = KK + [ks1 ks2 0 ;ks3 ks4 0 ;0 0 0];
    ks = kspring + gamma .* h ./ (beta .* h.^2) .* kc + 1 ./ (beta .* h.^2) .* MM; % Linearization

end
% end
