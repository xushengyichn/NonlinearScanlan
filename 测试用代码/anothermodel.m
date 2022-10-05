%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: Shengyi Xu xushengyichn@outlook.com
%Date: 2022-06-23 17:15:29
%LastEditors: Shengyi Xu xushengyichn@outlook.com
%LastEditTime: 2022-06-26 12:03:18
%FilePath: \twindeck_ID\test.m
%Description: 
%
%Copyright (c) 2022 by Shengyi Xu xushengyichn@outlook.com, All Rights Reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear
close all

a=0.5;
b=0.5;
c=0.1;
d=20;
e=5;
t=-100:0.1:100;
y=a*exp(-exp(b-c.*(t-d)))+e;
plot(t,y)



