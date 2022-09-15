%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%Date: 2022-09-14 10:24:48
%LastEditors: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%LastEditTime: 2022-09-14 10:33:50
%FilePath: \NonlinearScanlan\ReadExpData.m
%Description: 本函数实现输入实验名，读取实验数据，同时返回上下游梁的振动数据
%
%Copyright (c) 2022 by xushengyichn 54436848+xushengyichn@users.noreply.github.com, All Rights Reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [up_t, UP, down_t, DOWN] = ReadExpData(ExpName)
    fname = ExpName;
    chanum = 8; %记录通道数

    for i = 1:chanum
        refname = strcat(fname, '#', int2str(i), '.tsp');
        fid = fopen(refname, 'r');
        fsamp = fscanf(fid, '%f,%*d,%*d,%*d,%*d,%*f,%*d,%f'); %只读取tsp文件的最后一个参数
        up_amp(i) = fsamp(2, 1);
        fs = fsamp(1, 1);
        %amp(i)为DASP中的通道传感器的标定值  amp(i)
        fclose(fid);
    end

    for j = 1:chanum
        dataname = strcat(fname, '#', int2str(j), '.sts');
        fid = fopen(dataname, 'rb');
        [DATA(:, j), cn] = fread(fid, 'single'); %读取二进制文件，新买的DASP
        fclose(fid);
        DATA(:, j) = DATA(:, j) ./ 50; %进制转换，新买的dasp        10000mv， 10v对应量程20cm，200mm
        D(:, j) = DATA(:, j);
    end

    up_D = D(:, 1:4); %
    down_D = D(:, 5:8);
    n = length(D(:, 1));

    for i = 1:n
        AM(i, 1) = 1 / fs * i;
    end

    t = AM(:, 1); %未经过裁切的时间序列

    up_index_start = 1;
    up_index_end = length(t);
    down_index_start = 1;
    down_index_end = length(t);
    disp("Windward girder: Start point: " + num2str(up_index_start) + ";End point: " + num2str(up_index_end) + ";")
    disp("Leeward girder: Start point: " + num2str(down_index_start) + ";End point: " + num2str(down_index_end) + ";")
    up_D = up_D(up_index_start:up_index_end, :);
    down_D = down_D(down_index_start:down_index_end, :);
    up_t = t(up_index_start:up_index_end);
    down_t = t(down_index_start:down_index_end);
    up_t = up_t - up_t(1);
    down_t = down_t - down_t(1);

    UP(:, 1) = (up_D(:, 1) - mean(up_D(:, 1)) + up_D(:, 2) - mean(up_D(:, 2)) + up_D(:, 3) - mean(up_D(:, 3)) + up_D(:, 4) - mean(up_D(:, 4))) / 4/1000; %上游梁振动响应（m）
    DOWN(:, 1) = (down_D(:, 1) - mean(down_D(:, 1)) + down_D(:, 2) - mean(down_D(:, 2)) + down_D(:, 3) - mean(down_D(:, 3)) + down_D(:, 4) - mean(down_D(:, 4))) / 4/1000; %下游梁振动响应（m）

end
