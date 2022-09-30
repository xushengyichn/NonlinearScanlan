%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%Date: 2022-09-29 12:05:24
%LastEditors: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%LastEditTime: 2022-09-30 14:50:52
%FilePath: \NonlinearScanlan\polynomial_zeta.m
%Description: 生成多项式模型气动阻尼随振幅变化曲线
%
%Copyright (c) 2022 by xushengyichn 54436848+xushengyichn@users.noreply.github.com, All Rights Reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function    [Zeta]=polynomial_zeta(Amplitude,a1,a2,a3,a4,a5,rho,U,D,omega0,m,mode_integral_1,mode_integral_2,mode_integral_3,mode_integral_4,mode_integral_5)
    Zeta=-rho*U*D*(a1*mode_integral_1+4*a2.*Amplitude/3/pi*mode_integral_2+a3.*Amplitude.^2/4*mode_integral_3+8*a4.*Amplitude.^3/15/pi^2*mode_integral_4+a5.*Amplitude.^4/8*mode_integral_5)/2/omega0/m;
end