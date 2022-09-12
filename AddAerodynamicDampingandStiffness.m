%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%Date: 2022-09-12 11:31:13
%LastEditors: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%LastEditTime: 2022-09-12 11:48:10
%FilePath: \NonlinearScanlan\AddAerodynamicDampingandStiffness.m
%Description: 在节段模型矩阵中添加空气动力阻尼和气动刚度
%
%Copyright (c) 2022 by xushengyichn 54436848+xushengyichn@users.noreply.github.com, All Rights Reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [C,K] = AddAerodynamicDampingandStiffness(C,K,rho,U,D,a1,H4)
    C(1,1)=C(1,1)-rho*U*D*a1;  
    K(1,1)=K(1,1)-rho*U^2*H4;
end