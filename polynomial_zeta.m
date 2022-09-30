%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%Date: 2022-09-29 12:05:24
%LastEditors: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%LastEditTime: 2022-09-29 12:07:42
%FilePath: \NonlinearScanlan\polynomial_zeta.m
%Description: 生成多项式模型气动阻尼随振幅变化曲线
%
%Copyright (c) 2022 by xushengyichn 54436848+xushengyichn@users.noreply.github.com, All Rights Reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function    [Zeta]=polynomial_zeta(Amplitude,a1,a2,a3,a4,a5,rho,U,D,omega0,m)
    Zeta=-rho*U*D*(a1+4*a2.*Amplitude/3/pi+a3.*Amplitude.^2/4+8*a4.*Amplitude.^3/15/pi^2+a5.*Amplitude.^4/8)/2/omega0/m;
end