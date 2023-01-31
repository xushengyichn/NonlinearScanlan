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

function [result] = b_0_6_tmd(number_of_modes_to_control,number_of_modes_to_consider,number_of_tmds,modal_damping_ratios,t_length,mTMD1,mTMD2,mTMD3,mTMD4,mTMD5,fTMD1,fTMD2,fTMD3,fTMD4,fTMD5,fTMD6,dTMD1,dTMD2,dTMD3,dTMD4,dTMD5,dTMD6,xTMD1,xTMD2,xTMD3,xTMD4,xTMD5,xTMD6,total_tmd_mass)
%     number_of_modes_to_control=[1,2,3,4,5,6]
%     number_of_modes_to_consider=10
%     number_of_tmds=3
%     modal_damping_ratios=ones(1,number_of_modes_to_consider)*0.003
%     total_tmd_mass_ratio = 0.0002 % 总质量比 The total mass ratio
%     mass_six_span = 10007779.7 % 深中通道非通航桥六跨连续梁质量 The mass of 6-span continuous beam of the non-navigational bridge of the Zhenzhong-Link
%     total_tmd_mass = total_tmd_mass_ratio * mass_six_span % 总质量 The total mass
%     dTMD1=0.3134 ;
%     dTMD2= 0.6923;
%     dTMD3=0.8764;
%     fTMD1=0.8946;
%     fTMD2=0.08504;
%     fTMD3=0.03905;
%     mTMD1=0.1698;
%     mTMD2=0.8781;
%     xTMD1=0.09835;
%     xTMD2=0.4211 ;
%     xTMD3=0.957;
%     t_length=100;
%     mTMD1=0.01*total_tmd_mass+mTMD1*total_tmd_mass*0.49 % 质量 The mass mTMD1
%     mTMD2=0.01*total_tmd_mass+mTMD2*total_tmd_mass*0.49 % 质量 The mass mTMD2
%     mTMD3=total_tmd_mass-mTMD1-mTMD2 % 质量 The mass mTMD3
%     fTMD1=0.7+fTMD1*0.3 % 频率 The frequency fTMD1
%     fTMD2=0.7+fTMD2*0.3 % 频率 The frequency fTMD2
%     fTMD3=0.7+fTMD3*0.3 % 频率 The frequency fTMD3
%     dTMD1=0.02+dTMD1*0.18 % 阻尼比 The damping ratio dTMD1
%     dTMD2=0.02+dTMD2*0.18 % 阻尼比 The damping ratio dTMD2
%     dTMD3=0.02+dTMD3*0.18 % 阻尼比 The damping ratio dTMD3
%     xTMD1=xTMD1*660 % TMD1的x坐标 The x-coordinate of TMD1
%     xTMD2=xTMD2*660 % TMD2的x坐标 The x-coordinate of TMD2
%     xTMD3=xTMD3*660 % TMD3的x坐标 The x-coordinate of TMD3

    mTMD6=total_tmd_mass-mTMD1-mTMD2-mTMD3-mTMD4-mTMD5; % 质量 The mass mTMD6;
    TMDs_mass=[mTMD1,mTMD2,mTMD3,mTMD4,mTMD5,mTMD6];
    TMDs_frequency=[fTMD1,fTMD2,fTMD3,fTMD4,fTMD5,fTMD6];
    TMDs_damping_ratio=[dTMD1,dTMD2,dTMD3,dTMD4,dTMD5,dTMD6];
    TMDs_location=[xTMD1,xTMD2,xTMD3,xTMD4,xTMD5,xTMD6];
    [result] = a_0_main(number_of_modes_to_control,number_of_modes_to_consider,number_of_tmds,modal_damping_ratios,t_length,TMDs_mass,TMDs_frequency,TMDs_damping_ratio,TMDs_location);
    result=result.dis_all_modes_sum;

end
