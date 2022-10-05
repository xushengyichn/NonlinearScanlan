%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%Date: 2022-09-15 11:02:00
%LastEditors: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%LastEditTime: 2022-09-15 11:06:01
%FilePath: \NonlinearScanlan\CalData_Polynomial_noTMD_singledegree.m
%Description: 该函数用于计算多项式模型的振动响应
%
%Copyright (c) 2022 by xushengyichn 54436848+xushengyichn@users.noreply.github.com, All Rights Reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [t, dis] = CalData_Polynomial_noTMD_singledegree(ExpName, girderindex)
    %该函数用于计算多项式模型下安装TMD后的振动响应
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
    end

    % 结构参数
    % Structural parameters
    D = 0.667; % deck depth
    m = 80; % mass of the segment model
    Mass = m;
    L=3.6;%节段模型长度
    a=a*L;%识别的气动参数为单位长度下的，需要根据节段模型实际长度进行扩展（实际增加的是气动力，为了简化直接修正涡激力参数）。
    Fre= F0;
    rho = 1.225; % density of the air
    U = my_table.Windspeed(isexist); % density of the air
    

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

    %% Calculate the response
    u0 = [-0.086]*10e-3;
    udot0 = [0];
%     u0=-0.0023
%     udot0=0.0436
    P = zeros(1, length(tt));
    Lb=Buffeting(rho,t,0.02,U,0.078,0.015,-0.165);%见深中通道110m报告
    P(1,:) = Lb;
    nModes = 1;
    MM=m;
    CC=2*MM*(2*pi*Fre)*Zeta0;
    KK=MM*(2*pi*Fre)^2;

    out = polynomial_NB_adstiff_withlimit(Fre, Mass, Zeta0, rho, D, U, a,H4, t, P, u0, udot0,upperlimit,lowerlimit,Fre)

    dis=out(:,2);%节段模型位移
end