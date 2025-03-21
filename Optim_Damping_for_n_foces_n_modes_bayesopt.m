%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%Date: 2022-10-17 15:52:12
%LastEditors: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%LastEditTime: 2022-11-11 16:17:50
%FilePath: \NonlinearScanlan\Optim_Damping_for_n_foces_n_modes_bayesopt.m
%Description: 返回优化所需的多阶级模态最小阻尼
%
%Copyright (c) 2022 by xushengyichn 54436848+xushengyichn@users.noreply.github.com, All Rights Reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calmodes_all 考虑结构矩阵的n阶模态，如果n是固定的，输入一个数即可
function [minDamping_allmodes,result]=Optim_Damping_for_n_foces_n_modes_bayesopt(mode_numbers,numberofTMD,mTMD1,mTMD2,zetaTMD1,zetaTMD2,zetaTMD3,fTMD1,fTMD2,fTMD3,xTMD1,xTMD2,xTMD3,calmodes_all,mu)
    
    mass_six_span = 10007779.7;
    massallTMD=mass_six_span*mu;
    mTMD3=massallTMD-mTMD1-mTMD2;
    
    if mTMD3<=0
        penalty=1;
         mTMD3=1;
         mTMD1=1;
         mTMD2=1;
    else
        penalty=0;
    end


    mTMD=[mTMD1 mTMD2 mTMD3];
    zetaTMD=[zetaTMD1 zetaTMD2 zetaTMD3];
    fTMD=[fTMD1 fTMD2 fTMD3];
    xTMD=[xTMD1 xTMD2 xTMD3];
    mindampingRatio_modes=zeros(length(mode_numbers),1);
    mindampingRatio_modes_sys=zeros(length(mode_numbers),1);
    resultall_modes=[];
    for k1 = 1:length(mode_numbers) 
        mode_number=mode_numbers(k1);
        [minDampingRatio,minDampingRatio_sys,resultall]=Optim_Damping_for_n_modes(mode_number,numberofTMD,mTMD,zetaTMD,fTMD,xTMD,calmodes_all);
        mindampingRatio_modes(k1)=min(minDampingRatio);
        mindampingRatio_modes_sys(k1)=min(minDampingRatio_sys);
        str="resultall_modes.result"+num2str(k1)+"=resultall";
        eval(str)
    end
    result.mindampingRatio_modes=mindampingRatio_modes;
    result.mindampingRatio_modes_sys=mindampingRatio_modes_sys;
    result.resultall_modes=resultall_modes;

    if penalty==1
        minDamping_allmodes=min(mindampingRatio_modes);%设置一个非常大的数
    else
        minDamping_allmodes=min(mindampingRatio_modes);
    end

end