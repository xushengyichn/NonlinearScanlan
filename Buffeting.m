%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%Date: 2022-09-15 10:04:36
%LastEditors: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%LastEditTime: 2022-09-15 10:36:43
%FilePath: \NonlinearScanlan\Buffeting.m
%Description: 本函数用于生成抖振力升力时程，不考虑自激力《桥梁风工程》公式7-1-57~59
%
%Copyright (c) 2022 by xushengyichn 54436848+xushengyichn@users.noreply.github.com, All Rights Reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function    [Lb] = Buffeting(rho,t,I,U,Cl,Cl_deri,Cd)
%     t=0:1/256:60;
%     I=0.02;
%     U=5;
    sigma=I*U;
    u=normrnd(0,sigma,1,length(t));
    w=normrnd(0,sigma,1,length(t));
    Lb=1/2*rho*U^2*(Cl*(2*u/U)+(Cl_deri+Cd)*w/U);
end

