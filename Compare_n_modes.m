%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%Date: 2022-10-13 10:09:09
%LastEditors: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%LastEditTime: 2022-10-17 14:12:24
%FilePath: \NonlinearScanlan\Compare_n_modes.m
%Description: 计算一阶模态与多阶模态下的阻尼比（函数）
%
%Copyright (c) 2022 by xushengyichn 54436848+xushengyichn@users.noreply.github.com, All Rights Reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% 求解安装TMD后的响应
function [result]=Compare_n_modes(mode_number,numberofTMD,mTMD,zetaTMD,fTMD,locationTMD,calmodes)
% numberofTMD 所需要优化的TMD数量
% mTMD TMD质量
% zetaTMD TMD阻尼比
% fTMD TMD频率
% locationTMD TMD位置
% calmodes 计算模态数
% mode_number 气动力施加的模态

addpath("函数/")

% 读取模态数据
modeinfo=load('modeinfo.mat');
MM_eq=modeinfo.MM_eq(1:calmodes,1:calmodes);
KK_eq=modeinfo.KK_eq(1:calmodes,1:calmodes);
eig_val=modeinfo.eig_val(1:calmodes,1:calmodes);
eig_vec=modeinfo.eig_vec(:,1:calmodes);

%% 计算不安装TMD情况下各阶模态各点最大位移
% nTMD = 0;
% mTMD = 0;
% zetaTMD = 0;
% omegaTMD = 0;
% nodeTMD = 0;
% mode_numbers = 1:1:1;
% ifcalmode = 3;
% 
% for mode_number = 1:length(mode_numbers)
%     [modemaxdis_single_noTMD(mode_number), usinglemax_noTMD(:, mode_number), uallmax_noTMD(:, mode_number)] = CalData_Polynomial_withTMD_multidegree(nTMD, mTMD, zetaTMD, omegaTMD, nodeTMD, mode_number, ifcalmode, MM_eq, KK_eq, calmodes, eig_val, eig_vec);
% end

% clear nTMD mTMD zetaTMD omegaTMD nodeTMD mode_numbers ifcalmode
%%
% 设置TMD参数

% mass_six_span = 10007779.7;
% mu = mu;
% mTMDall = mass_six_span * mu;

nTMD = numberofTMD;
mTMD=mTMD;
Ftmd = fTMD;
zetaTMD = zetaTMD;
omegaTMD = 2 * pi * Ftmd;
% cTMD = [2 * mTMD(1) * omegatmds * 0.05];
% disp(cTMD)


% mode_number=1;%气动力施加在第一阶模态
% mode_numbers = 1:1:1;
ifcalmode = 3;
nodeTMD=locationTMD;
% [modemaxdis_single,usinglemax,uallmax,u1,output]=CalData_Polynomial_withTMD_multidegree(nTMD,mTMD,zetaTMD,omegaTMD,nodeTMD,mode_number,ifcalmode,MM_eq,KK_eq,calmodes,eig_val,eig_vec,t_length);
[result]=CalDamping_Polynomial_withTMD_multidegree(nTMD,mTMD,zetaTMD,omegaTMD,nodeTMD,mode_number,ifcalmode,MM_eq,KK_eq,calmodes,eig_val,eig_vec);
Mode_damping=result.Mode.("Damping ratio");
% disp(result.Mode)
% disp(result.Mode_sys)
% phi_nodeTMD=modes(find(modes(:,2)==nodeTMD),3:2+calmodes);%安装TMD位置的振型向量大小
% 
% for k1 = 1:calmodes
%     dis_nodeTMD(k1,:)=u1(k1,:).*phi_nodeTMD(k1);%计算每一阶模态在安装TMD位置的位移
% end
% % dis_nodeTMD_sum=sum(dis_nodeTMD);%计算安装TMD位置的位移
% dis_TMD=u1(end,:);%计算TMD的位移
% result.dis_TMD=dis_TMD;
% result.dis_nodeTMD=dis_nodeTMD;
% result.output=output;

end