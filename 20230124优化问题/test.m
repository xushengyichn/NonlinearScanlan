%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%Date: 2023-01-25 10:39:30
%LastEditors: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%LastEditTime: 2023-01-25 15:28:37
%FilePath: \20230124优化问题\a_0_main.m
%Description: 本函数为主函数，分别调用各个子函数完成桥梁在多阶级涡振下的响应总和计算
%
%Copyright (c) 2023 by xushengyichn 54436848+xushengyichn@users.noreply.github.com, All Rights Reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [result] = test(number_of_modes_to_control,number_of_modes_to_consider,number_of_tmds,modal_damping_ratios,t_length,mTMD1,mTMD2,fTMD1,fTMD2,fTMD3,dTMD1,dTMD2,dTMD3,xTMD1,xTMD2,xTMD3,total_tmd_mass)
    [result] = sum(mTMD1+mTMD2);
end
