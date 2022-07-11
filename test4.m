%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: Shengyi Xu xushengyichn@outlook.com
%Date: 2022-06-27 18:40:52
%LastEditors: Shengyi Xu xushengyichn@outlook.com
%LastEditTime: 2022-06-27 19:16:22
%FilePath: \NonlinearScanlan\polynomial_NB_withTMD.m
%Description: 本函数目的为计算多项式模型安装tmd后的的响应。
%
%Copyright (c) 2022 by Shengyi Xu xushengyichn@outlook.com, All Rights Reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This is a function
function out = test4(Fre, Mass, Zeta0, rho, D, U, a, t,h, P,  u0, udot0,nModes,matrixsize,mtmd)

    % Nonlinear Newmark's Direct Integration Method with polynomial model
    % (n = number of time steps)
    % (ndof = number degrees of freedom)

    % INPUT
    % Fre      = Frequency of the system         =>[1]
    % Mass     = Mass of the system              =>[1]
    % Zeta0    = Damping ratio of the system     =>[1]
    % rho      = Air density                     =>[1]
    % D        = Reference length                =>[1]
    % U        = Wind speed at a certain reduced frequency       =>[1]
    % a        = 多项式模型系数（四次多项式5个系数）       =>[1,5]
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
