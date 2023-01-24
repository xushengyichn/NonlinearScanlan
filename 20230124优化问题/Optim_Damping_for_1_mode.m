%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%Date: 2022-10-17 15:52:12
%LastEditors: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%LastEditTime: 2022-10-17 16:34:15
%FilePath: \NonlinearScanlan\Optim_Damping_for_n_modes.m
%Description: 返回优化所需的多阶级模态最小阻尼
%
%Copyright (c) 2022 by xushengyichn 54436848+xushengyichn@users.noreply.github.com, All Rights Reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [minDampingRatio,resultall]=Optim_Damping_for_1_mode(numberofTMD,mTMD,zetaTMD,fTMD,xTMD,calmodes_all)
% numberofTMD  TMD数量  N
% mTMD TMD质量 1*N
% zetaTMD TMD阻尼比 1*N
% fTMD TMD频率 1*N
% xTMD TMD位置 1*N
% calmodes_all 计算的气动力模态 1*M eg.[1 3] 这就表示计算第一阶和第二阶模态 
    minDamp_aero=zeros(length(calmodes_all),1);
    resultall=[];
    for k1 = 1: length(calmodes_all)
        calmode=calmodes_all(k1);
        [min_Mode_damping,result]=Compare_1_mode(numberofTMD,mTMD,zetaTMD,fTMD,xTMD,calmode);
       
        minDamp_aero(k1)=min_Mode_damping;
        str="resultall.result"+num2str(k1)+"=result";
        eval (str)
    end
    resultall.minDamp_aero=minDamp_aero;
    if sum(mTMD)>10007779.7*0.015*3
        minDampingRatio=-1;
    else
        minDampingRatio=min(minDamp_aero);
    end

end