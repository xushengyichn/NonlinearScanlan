clc;clear;close all

% 打开文件
% fid = fopen('logs.json', 'r');

fid = fopen('logs.json', 'r');
tline = fgetl(fid);
target =[];
while ischar(tline)
    data = jsondecode(tline);
    % 处理数据
    target=[target data.target];
    tline = fgetl(fid);
end
fclose(fid);
target(1)=[];
max(target)