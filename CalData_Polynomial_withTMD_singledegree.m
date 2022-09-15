%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%Date: 2022-09-14 10:49:00
%LastEditors: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%LastEditTime: 2022-09-14 11:23:29
%FilePath: \NonlinearScanlan\CalData_Polynomial_withTMD_singledegree.m
%Description: 该函数用于计算多项式模型下安装TMD后的振动响应
%
%Copyright (c) 2022 by xushengyichn 54436848+xushengyichn@users.noreply.github.com, All Rights Reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [t, dis] = CalData_Polynomial_withTMD_singledegree(ExpName, girderindex, TMDsindex, freshift)
    %该函数用于计算多项式模型下安装TMD后的振动响应
    %ExpName:实验名称
    %girderindex:梁号(上游输入1，下游输入2)
    %TMDsindex:安装TMD的需要（TMD的序号，每个TMD对应的序号见TMD_logfile.mat）
    %freshift:频率偏移，因为TMD频率分辨率很低，在一定范围中移动可以更好与试验结果相对照

    % 导入试验数据
    my_table_tmd=load('TMD_logfile.mat');
    my_table_tmd=my_table_tmd.my_table;
    my_table=load('SZTD110_logfile.mat');
    my_table=my_table.my_table;

%     ExpName='SZTD-110-case2-22.3-fasan-2401';
%     girderindex=1;
%     TMDsindex=[12 14 15 18];
%     freshift=0.16;

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
    
    end

    % 结构参数
    % Structural parameters
    D = 0.667; % deck depth
    m = 80; % mass of the segment model
    Mass = m;
    F0 = my_table.up_Fre_vibration(isexist); % Frequency without wind
    Fre= F0;
    rho = 1.225; % density of the air
    U = my_table.Windspeed(isexist); % density of the air
    Zeta0 =  my_table.up_dltx_zeta0(isexist); % damping ratio without wind

    h = 1/256; % 时间步长 % Step size of the algorithm
    t=0:h:60; % Time vector
    told = t;
    tt = told;
    P = zeros(1, size(told, 1));

    disp("节段模型质量："+num2str(Mass));
    disp("节段模型阻尼系数："+num2str(Zeta0));
    disp("节段模型频率："+num2str(Fre));
    disp("节段模型风速："+num2str(U));
    disp("节段模型气动阻尼参数"+num2str(a));
    disp("节段模型气动刚度参数"+num2str(H4));

    %% TMD 参数
    sel=TMDsindex;  
    zetatmd = my_table_tmd.zeta(sel);
    fretmd = my_table_tmd.fre(sel)+freshift;
    mtmd = ones(length(sel),1)*0.25;
    disp("TMD质量："+num2str(mtmd));
    disp("TMD阻尼系数："+num2str(zetatmd));
    disp("TMD频率："+num2str(fretmd));

    %% Calculate the response
    u0 = [-0.086; -1;-1;-1;-1]*10e-3;
    udot0 = [0; 0;0;0;0];

    P = zeros(5, length(tt));
    nModes = 1;
    matrixsize=length(TMDsindex);

    [MM,CC,KK]=CreateMatrixwithTMD(nModes,Mass,Zeta0,Fre,mtmd,zetatmd,fretmd);

    a1_lower=a1+4/3*a2/pi*lowerlimit+a3/4*lowerlimit^2+8/15*a4/pi*lowerlimit^3+a5/8*lowerlimit^4;
    disp("lowerlimit="+a1_lower)

    [CC1,KK1]=AddAerodynamicDampingandStiffness(CC,KK,rho,U,D,a1_lower,H4);
    [V,DD]=eigs(KK1,MM);
    Result0=sort(diag(sqrt(DD)/2/pi));
    Result1= Complex_Eigenvalue_Analysis(MM,CC,KK1);%不考虑气动阻尼的复数特征值分析
    Result2=Complex_Eigenvalue_Analysis(MM,CC1,KK1);%考虑气动阻尼的复特征值分析
    disp(Result2)
    out = test_polynomial_NB_withTMDs_addstiff_withlimit(Fre, Mass, Zeta0, rho, D, U, a,H4, t,h, P, u0, udot0,upperlimit,lowerlimit,Fren_vibration_withwind,nModes,mtmd,fretmd,zetatmd);
    dis=out(:,2);%节段模型位移
end