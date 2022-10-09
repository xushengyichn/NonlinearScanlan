%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%Date: 2022-10-09 17:47:03
%LastEditors: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%LastEditTime: 2022-10-09 17:50:37
%FilePath: \NonlinearScanlan\CalData_Polynomial_withTMD_multidegree_multifocemode.m
%Description: 计算多阶气动力施加后，每个节点最大位移的和
%
%Copyright (c) 2022 by xushengyichn 54436848+xushengyichn@users.noreply.github.com, All Rights Reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [modemaxdis_single,usinglemax,uallmax,totalmax]=CalData_Polynomial_withTMD_multidegree_multifocemode(nTMD,mTMD,zetaTMD,omegaTMD,nodeTMD,mode_numbers,ifcalmode,MM_eq,KK_eq,calmodes,eig_val,eig_vec)
% function [uallmax]=CalData_Polynomial_withTMD_multidegree_multifocemode(nTMD,mTMD,zetaTMD,omegaTMD,nodeTMD,mode_numbers,ifcalmode,MM_eq,KK_eq,calmodes,eig_val,eig_vec)
    for k1 = 1:mode_numbers
            %计算每个模态的最大位移
            [modemaxdis_single(k1),usinglemax(k1,:),uallmax(k1,:)] = CalData_Polynomial_withTMD_multidegree(nTMD,mTMD,zetaTMD,omegaTMD,nodeTMD,k1,ifcalmode,MM_eq,KK_eq,calmodes,eig_val,eig_vec);
    end
    totalmax=sum(sum(uallmax));%所有模态所有节点位移之和
end