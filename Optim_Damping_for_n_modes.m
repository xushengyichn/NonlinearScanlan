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


function [minDampingRatio,minDampingRatio_sys,resultall]=Optim_Damping_for_n_modes(mode_number,numberofTMD,mTMD,zetaTMD,fTMD,xTMD,calmodes_all)
    minDamp_sys=zeros(length(calmodes_all),1);
    minDamp_aero=zeros(length(calmodes_all),1);
    resultall=[];
    for k1 = 1: length(calmodes_all)
        calmodes=calmodes_all(k1);
        [result]=Compare_n_modes(mode_number,numberofTMD,mTMD,zetaTMD,fTMD,xTMD,calmodes);
        minDamp_sys(k1)=min(result.Mode_sys.("Damping ratio"));
        minDamp_aero(k1)=min(result.Mode.("Damping ratio"));
        str="resultall.result"+num2str(k1)+"=result";
        eval(str)
    end
    minDampingRatio=min(minDamp_aero);
    minDampingRatio_sys=min(minDamp_sys);
end