%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%Date: 2022-10-13 11:11:22
%LastEditors: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%LastEditTime: 2022-10-13 11:12:36
%FilePath: \NonlinearScanlan\20221012考虑一阶模态与多阶区别\Main.m
%Description: 主文件（调用对比函数进行批量分析）
%
%Copyright (c) 2022 by xushengyichn 54436848+xushengyichn@users.noreply.github.com, All Rights Reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc;clear;close all;

% mu=0.05/100;
% zetaTMD=0.06;
% fTMD=0.8339;
% locationTMD=1036;
% calmodes=2;
% t_length=200;
% [result]=Compare_n_modes(mu,zetaTMD,fTMD,locationTMD,calmodes,t_length);
% disnode=result.dis_nodeTMD;
% disnode_1=disnode(1,:);
% disnode_2=disnode(2,:);
% disnode_sum=disnode_1+disnode_2;
% figure
% fs=100;
% t=0:1/fs:t_length;
% plot(t,disnode_1)
% hold on 
% plot(t,disnode_2)
% plot(t,disnode_sum)
% legend('1','2','sum')
% data=disnode_sum';
% [psd_avg, f, psd_plot] = fft_transfer(fs,disnode_sum');
% figure
% plot(f,psd_plot)
% [psd_avg, f, psd_plot] = fft_transfer(fs,disnode_1');
% figure
% plot(f,psd_plot)
% [psd_avg, f, psd_plot] = fft_transfer(fs,disnode_2');
% figure
% plot(f,psd_plot)


mus=0.01/100:0.01/100:5/100;
Frequency=zeros(calmodes+1,length(mus));
Damping_ratio=zeros(calmodes+1,length(mus));
for k1 = 1:length(mus)
    mu=mus(k1);
    zetaTMD=0.06;
    fTMD=0.8339;
    locationTMD=1036;
    calmodes=2;
    t_length=200;
    [result]=Compare_n_modes(mu,zetaTMD,fTMD,locationTMD,calmodes,t_length);
    disnode=result.dis_nodeTMD;
    disnode_1=disnode(1,:);
    disnode_2=disnode(2,:);
    disnode_sum=disnode_1+disnode_2;
    Frequency(:,k1)=result.output.Mode(:,2);
    Damping_ratio(:,k1)=result.output.Mode(:,3);
end
% disnode=result.dis_nodeTMD;
% disnode_1=disnode(1,:);
% disnode_2=disnode(2,:);
% disnode_sum=disnode_1+disnode_2;
% figure
% fs=100;
% t=0:1/fs:t_length;
% plot(t,disnode_1)
% hold on 
% plot(t,disnode_2)
% plot(t,disnode_sum)
% legend('1','2','sum')
% data=disnode_sum';
% [psd_avg, f, psd_plot] = fft_transfer(fs,disnode_sum');
% figure
% plot(f,psd_plot)
% [psd_avg, f, psd_plot] = fft_transfer(fs,disnode_1');
% figure
% plot(f,psd_plot)
% [psd_avg, f, psd_plot] = fft_transfer(fs,disnode_2');
% figure
% plot(f,psd_plot)