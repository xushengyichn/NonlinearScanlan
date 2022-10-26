%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%Date: 2022-10-14 11:41:51
%LastEditors: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%LastEditTime: 2022-10-26 17:57:22
%FilePath: \NonlinearScanlan\Optimization_Damping.m
%Description: 直接计算安装TMD后的阻尼比
%
%Copyright (c) 2022 by xushengyichn 54436848+xushengyichn@users.noreply.github.com, All Rights Reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc;clear;close all;
addpath("函数\")
calmodes=1;

% 
% [Ftmds_M,zetaTMDs_M,Massratio_M]=ndgrid(Frequency,Damping_ratio,mus);
% variables = [Ftmds_M(:),zetaTMDs_M(:),Massratio_M(:)];
% 
% 
% numIterations = size(variables,1);
% res = zeros(numIterations, 1);


%% 单点计算
if 0
mass_six_span = 10007779.7;
mode_number=1;%气动力施加在第一阶模态

% numberofTMD=1;
% mTMD=[mass_six_span*0.015];
% zetaTMD=[0.1];
% fTMD=[0.82];
% xTMD=56.2
load opt_1tmd_2modes
mTMD=sol.mTMD;
zetaTMD=sol.zetaTMD;
fTMD=sol.fTMD;
xTMD= sol.xTMD;


calmodes=2;
calmodes_all=[1];


[result]=Compare_n_modes(mode_number,numberofTMD,mTMD,zetaTMD,fTMD,xTMD,calmodes);
result.Mode
result.Mode_sys
% tic
% [minDampingRatio,minDampingRatio_sys]=Optim_Damping_for_n_modes(mode_number,numberofTMD,mTMD,zetaTMD,fTMD,xTMD,calmodes_all);
% toc
% disp(minDampingRatio)
% disp(minDampingRatio_sys)
end

%% 穷举问题(阻尼、频率、质量比)
if 0
mode_number=1;%气动力施加在第一阶模态
numberofTMD=1;
mass_six_span = 10007779.7;
zetaTMDs=0.05:0.001:0.25;
fTMDs=0.7:0.001:0.9;
muTMDs=0.003:0.003:0.03;
mTMDs=muTMDs*mass_six_span;
xTMDs=56.2;

[zetaTMDs_M,fTMDs_M,xTMDs_M]=ndgrid(zetaTMDs,fTMDs,mTMDs);
variables = [zetaTMDs_M(:),fTMDs_M(:),xTMDs_M(:)];
numIterations=size(variables,1);
ppm = ParforProgressbar(size(variables,1));
parfor k1 = 1:numIterations
% for k1 = 1:numIterations
    zetaTMD=variables(k1,1);
    fTMD=variables(k1,2);
    mTMD=variables(k1,3);
    xTMD=xTMDs;
    [minDampingRatio]=Optim_Damping_for_n_modes(mode_number,numberofTMD,mTMD,zetaTMD,fTMD,xTMD,calmodes);
    res(k1,1)=minDampingRatio;
%     disp(res)
    pause(100/numIterations);
    ppm.increment();
end

result=[variables res];

str="save traverse_"+num2str(numberofTMD)+"TMD_"+num2str(calmodes)+"modes_xloc"+num2str(xTMDs)+".mat result";
eval(str)

end
%% 数据分析
load traverse_1TMD_1modes_xloc56.2

Damping_result_filter_data=result;
Damping_result_filter_data=result(result(:,4)>=0,:);
m0=10007779.7;
f0=0.833853594612215;
Damping_result_filter_data(:,5)=Damping_result_filter_data(:,2)/f0;
Damping_result_filter_data(:,6)=Damping_result_filter_data(:,3)/m0;

a=zeros(size(Damping_result_filter_data,1),1);
b=zeros(size(Damping_result_filter_data,1),1);
c=zeros(size(Damping_result_filter_data,1),1);

Damping_result_filter_data(:,6)=round(Damping_result_filter_data(:,6),3);
fixmu=Damping_result_filter_data(Damping_result_filter_data(:,6)==0.015,:);

% for k1 =1:size(Damping_result_filter_data,1)
%      if mod(Damping_result_filter_data(k1,1)*1000,2)==0
%          a(k1)=1;
%      end
%      if mod(Damping_result_filter_data(k1,2)*1000,2)==0
%          b(k1)=1;
%      end
%      if and(a(k1)==1, b(k1)==1)
%          c(k1)=1;
%      end
% end
% Damping_result_filter_data(c==0,:)=[];


figure
scatter(fixmu(:,5),fixmu(:,1))
figure
scatter3(Damping_result_filter_data(:,1),Damping_result_filter_data(:,5),Damping_result_filter_data(:,6),[],Damping_result_filter_data(:,4))

save Dampingresult_filter_data_traverse_1TMD_1modes_xloc56.2.mat Damping_result_filter_data fixmu

%% 穷举问题(阻尼、频率、位置)
if 0
mode_number=1;%气动力施加在第一阶模态
numberofTMD=1;
mass_six_span = 10007779.7;
zetaTMDs=0.05:0.001:0.25;
fTMDs=0.8:0.001:1.0;
mTMDs=0.015*mass_six_span;
xTMDs=0:1:606;
[zetaTMDs_M,fTMDs_M,xTMDs_M]=ndgrid(zetaTMDs,fTMDs,xTMDs);
variables = [zetaTMDs_M(:),fTMDs_M(:),xTMDs_M(:)];
numIterations=size(variables,1);
% ppm = ParforProgressbar(size(variables,1));
% parfor k1 = 1:numIterations
for k1 = 1:numIterations
    zetaTMD=variables(k1,1);
    fTMD=variables(k1,2);
    xTMD=variables(k1,3);
    [minDampingRatio]=Optim_Damping_for_n_modes(mode_number,numberofTMD,mTMDs,zetaTMD,fTMD,xTMD,calmodes);
    res(k1)=minDampingRatio;
%     disp(res)
%     pause(100/numIterations);
%     ppm.increment();
end
end



%% 优化问题(单模态一个TMD)
if 0
mode_number=1;%气动力施加在第一阶模态
numberofTMD=1;
mass_six_span = 10007779.7;
prob = optimproblem('Description','Optimize the TMDs','ObjectiveSense','maximize'); %优化为最优阻尼比
zetaTMD =optimvar('zetaTMD',numberofTMD,'LowerBound',0.05,'UpperBound',0.25); %阻尼比
fTMD =optimvar('fTMD',numberofTMD,'LowerBound',0.8,'UpperBound',1.0); %频率
mTMD =optimvar('mTMD',numberofTMD,'LowerBound',0.015*mass_six_span,'UpperBound',0.015*mass_six_span); %质量
xTMD = optimvar('xTMD', numberofTMD,'LowerBound', 56.2, 'UpperBound', 56.2); %TMD所安装的节点
calmodes_all=1:1:1;
% options = optimoptions('particleswarm',Display='iter',PlotFcn='pswplotbestf');
options = optimoptions('particleswarm',Display='iter',PlotFcn='pswplotbestf',UseParallel=true);
% options = optimoptions('ga', 'Display', 'iter', 'PlotFcn', {'gaplotscorediversity', 'gaplotbestf', 'gaplotrange'}, 'UseParallel', true);
[minDampingRatio,minDampingRatio_sys] = fcn2optimexpr(@Optim_Damping_for_n_modes,mode_number,numberofTMD,mTMD,zetaTMD,fTMD,xTMD,calmodes_all);

% msum = sum(mTMD) == mass_six_span*0.015;
% prob.Constraints.msum = msum;

% prob.Constraints.fTMD = fTMD ==0.82;
% prob.Constraints.locationTMD = locationTMD== 1179;

prob.Objective=minDampingRatio;
show(prob)

% x0.zetaTMD = 0.1*ones(numberofTMD,1);
% x0.fTMD = 0.82*ones(numberofTMD,1);
% x0.mTMD = 0.015*mass_six_span*ones(numberofTMD,1);
% x0.xTMD = 0.1*ones(numberofTMD,1);

[sol, optval] = solve(prob,'Solver','particleswarm','Options',options);

val = evaluate(prob.Objective,sol);
disp(val)

strmdoes=(num2str(calmodes_all));
strmdoes=strjoin(strsplit(strmdoes),'_');

str="save opt_"+num2str(numberofTMD)+"TMD_"+strmdoes+"modes_xloc"+num2str(sol.xTMD)+".mat";
eval(str)
% save opt_1tmd_1modes
end
%% 数据分析
if 0
str="load opt_"+num2str(numberofTMD)+"TMD_"+strmdoes+"modes_xloc"+num2str(sol.xTMD)+".mat";
eval(str)

mTMD=sol.mTMD;
zetaTMD=sol.zetaTMD;
fTMD=sol.fTMD;
xTMD= sol.xTMD;
[result]=Optim_Damping_for_n_modes(mode_number,numberofTMD,mTMD,zetaTMD,fTMD,xTMD,calmodes_all);
result
calmodes=2
[result]=Compare_n_modes(mode_number,numberofTMD,mTMD,zetaTMD,fTMD,xTMD,calmodes);
result.Mode
result.Mode_sys
end


%% 两个TMD
if 0
a=importdata("opt_1TMD_1_2modes_xloc56.2.mat")
b=importdata("opt_1TMD_1modes_xloc56.2.mat")
mode_number=1
numberofTMD=2
mTMD=[a.sol.mTMD b.sol.mTMD];
zetaTMD=[a.sol.zetaTMD b.sol.zetaTMD];
fTMD=[a.sol.fTMD b.sol.fTMD];
xTMD=[56.2 289.465]
calmodes=2
[result]=Compare_n_modes(mode_number,numberofTMD,mTMD,zetaTMD,fTMD,xTMD,calmodes);
result.Mode
result.Mode_sys


% test.fTMD=0.93
% test.mTMD=sol.mTMD;
% test.xTMD=56.2
% test.zetaTMD=0.22
% evaluate(prob.Objective,test);


% load damping_1tmd_5modes.mat
% 
% [~,seq]=max(result(:,4));
% opt_design=result(seq,:)
end

%% 单模态一个1个TMD（基于1 2 3 阶气动力）
calmodes_alls=[1 2 3];
mode_numbers=[1 2 3];

result=[];
for k1 = 1:3
numberofTMD=1;
mass_six_span = 10007779.7;
prob = optimproblem('Description','Optimize the TMDs','ObjectiveSense','maximize'); %优化为最优阻尼比
zetaTMD =optimvar('zetaTMD',numberofTMD,'LowerBound',0.04,'UpperBound',0.25); %阻尼比
fTMD =optimvar('fTMD',numberofTMD,'LowerBound',0.6,'UpperBound',1.5); %频率
mTMD =optimvar('mTMD',numberofTMD,'LowerBound',0.015*mass_six_span,'UpperBound',0.015*mass_six_span); %质量
xTMD = optimvar('xTMD', numberofTMD,'LowerBound', 0, 'UpperBound', 660); %TMD所安装的节点
calmodes_all=calmodes_alls(k1);
options = optimoptions('particleswarm',Display='iter',PlotFcn='pswplotbestf',UseParallel=true);
[min_Mode_damping] = fcn2optimexpr(@Compare_1_mode,numberofTMD,mTMD,zetaTMD,fTMD,xTMD,calmodes_all);


prob.Objective=min_Mode_damping;
show(prob)


[sol, optval] = solve(prob,'Solver','particleswarm','Options',options);

val = evaluate(prob.Objective,sol);
disp(val)

strmdoes=(num2str(calmodes_all));
strmdoes=strjoin(strsplit(strmdoes),'_');

str="result.sol"+num2str(k1)+"=sol";
eval(str)
str="result.optval"+num2str(k1)+"=optval";
eval(str)
str="result.prob"+num2str(k1)+"=prob";
eval(str)
end
save optim_1tmd_123mode_separate result
%% 数据分析
load optim_1tmd_123mode_separate result
load modeinfo.mat
% 模态1
xTMD=result.sol1.xTMD;
figure
plot(modes(:,1),modes(:,3))
hold on 

[~,index]=sort(abs(nodegap-xTMD));%查找与xTMD最接近的点的排序
xResult=nodegap(index(1:2));%获取最接近的两个点的x坐标
mode2nodes=mode(index(1:2),1:1);%获取两个点坐标的y值
phi_result=interp1(xResult,mode2nodes,xTMD,'linear','extrap');%插值以后任意点的振型
scatter(xTMD,phi_result)
% 模态2
xTMD=result.sol2.xTMD;
figure
plot(modes(:,1),modes(:,4))
hold on 

[~,index]=sort(abs(nodegap-xTMD));%查找与xTMD最接近的点的排序
xResult=nodegap(index(1:2));%获取最接近的两个点的x坐标
mode2nodes=mode(index(1:2),2);%获取两个点坐标的y值
phi_result=interp1(xResult,mode2nodes,xTMD,'linear','extrap');%插值以后任意点的振型
scatter(xTMD,phi_result)

% 模态3

xTMD=result.sol3.xTMD;
figure
plot(modes(:,1),modes(:,5))
hold on 

[~,index]=sort(abs(nodegap-xTMD));%查找与xTMD最接近的点的排序
xResult=nodegap(index(1:2));%获取最接近的两个点的x坐标
mode2nodes=mode(index(1:2),3);%获取两个点坐标的y值
phi_result=interp1(xResult,mode2nodes,xTMD,'linear','extrap');%插值以后任意点的振型
scatter(xTMD,phi_result)



%% 单模态一个3个TMD（基于1 2 3 阶气动力）


calmodes_alls=[1 2 3]; %以单模态的形式计算3阶涡振控制


result=[];

numberofTMD=3;
mass_six_span = 10007779.7;
prob = optimproblem('Description','Optimize the TMDs','ObjectiveSense','maximize'); %优化为最优阻尼比
zetaTMD =optimvar('zetaTMD',numberofTMD,'LowerBound',0.04,'UpperBound',0.25); %阻尼比
fTMD =optimvar('fTMD',numberofTMD,'LowerBound',0.6,'UpperBound',1.5); %频率
mTMD =optimvar('mTMD',numberofTMD,'LowerBound',0.001*mass_six_span,'UpperBound',0.015*mass_six_span); %质量
xTMD = optimvar('xTMD', numberofTMD,'LowerBound', 0, 'UpperBound', 660); %TMD所安装的节点
calmodes_all=calmodes_alls;
options = optimoptions('particleswarm',Display='iter',PlotFcn='pswplotbestf',UseParallel=true);
[minDampingRatio] = fcn2optimexpr(@Optim_Damping_for_1_mode,numberofTMD,mTMD,zetaTMD,fTMD,xTMD,calmodes_all);


prob.Objective=minDampingRatio;
show(prob)


[sol, optval] = solve(prob,'Solver','particleswarm','Options',options);

val = evaluate(prob.Objective,sol);
disp(val)

strmdoes=(num2str(calmodes_all));
strmdoes=strjoin(strsplit(strmdoes),'_');
result.sol=sol;
result.optval=optval;
result.prob=prob;

save optim_3tmd_123mode_all result

%% 数据分析
load optim_3tmd_123mode_all
load modeinfo.mat
% 模态1
xTMD=result.sol.xTMD;
figure
plot(modes(:,1),modes(:,3))
hold on
plot(modes(:,1),modes(:,4))
plot(modes(:,1),modes(:,5))
hold on 
for k1 = 1:3
[~,index]=sort(abs(nodegap-xTMD(k1)));%查找与xTMD最接近的点的排序
xResult=nodegap(index(1:2));%获取最接近的两个点的x坐标
mode2nodes=mode(index(1:2),k1);%获取两个点坐标的y值
phi_result(k1)=interp1(xResult,mode2nodes,xTMD(k1),'linear','extrap');%插值以后任意点的振型
end
scatter(xTMD,phi_result)

tab=struct2table(result.sol);
disp(tab)

%% 数据分析对比单模态，一阶一阶设计或三阶一起设计
load modeinfo.mat
allresult=importdata("optim_3tmd_123mode_all.mat")
sepresult=importdata("optim_1tmd_123mode_separate.mat")

figure
plot(modes(:,1),modes(:,3))
hold on 
plot(modes(:,1),modes(:,4))
plot(modes(:,1),modes(:,5))

% 模态1
xTMD=sepresult.sol1.xTMD;
[~,index]=sort(abs(nodegap-xTMD));%查找与xTMD最接近的点的排序
xResult=nodegap(index(1:2));%获取最接近的两个点的x坐标
mode2nodes=mode(index(1:2),1:1);%获取两个点坐标的y值
phi_result=interp1(xResult,mode2nodes,xTMD,'linear','extrap');%插值以后任意点的振型
scatter(xTMD,phi_result)
phi_result1=phi_result;
% 模态2
xTMD=sepresult.sol2.xTMD;
[~,index]=sort(abs(nodegap-xTMD));%查找与xTMD最接近的点的排序
xResult=nodegap(index(1:2));%获取最接近的两个点的x坐标
mode2nodes=mode(index(1:2),2);%获取两个点坐标的y值
phi_result=interp1(xResult,mode2nodes,xTMD,'linear','extrap');%插值以后任意点的振型
scatter(xTMD,phi_result)
phi_result2=phi_result;
% 模态3
xTMD=sepresult.sol3.xTMD;
[~,index]=sort(abs(nodegap-xTMD));%查找与xTMD最接近的点的排序
xResult=nodegap(index(1:2));%获取最接近的两个点的x坐标
mode2nodes=mode(index(1:2),3);%获取两个点坐标的y值
phi_result=interp1(xResult,mode2nodes,xTMD,'linear','extrap');%插值以后任意点的振型
scatter(xTMD,phi_result)
phi_result3=phi_result;

% 三个模态一起设计
fTMDall=allresult.sol.fTMD;
[r1,s1]=sort(fTMDall);%获取按照频率排序的序列
xTMDall=allresult.sol.xTMD;
xTMDall=xTMDall(s1);%按照频率顺序排列TMD为位置，为每个频率寻找对应的振型
for k1 = 1:3
[~,index]=sort(abs(nodegap-xTMDall(k1)));%查找与xTMD最接近的点的排序
xResult=nodegap(index(1:2));%获取最接近的两个点的x坐标
mode2nodes=mode(index(1:2),k1);%获取两个点坐标的y值
phi_result(k1)=interp1(xResult,mode2nodes,xTMDall(k1),'linear','extrap');%插值以后任意点的振型
end
scatter(xTMDall,phi_result,'filled')

sepresultcombine_fTMD=[sepresult.sol1.fTMD sepresult.sol2.fTMD sepresult.sol3.fTMD];
sepresultcombine_mTMD=[sepresult.sol1.mTMD sepresult.sol2.mTMD sepresult.sol3.mTMD];
sepresultcombine_xTMD=[sepresult.sol1.xTMD sepresult.sol2.xTMD sepresult.sol3.xTMD];
sepresultcombine_zetaTMD=[sepresult.sol1.zetaTMD sepresult.sol2.zetaTMD sepresult.sol3.zetaTMD];
sepresultcombine_phi_result=[phi_result1 phi_result2 phi_result3];
sepresultcombine.fTMD=sepresultcombine_fTMD;
sepresultcombine.mTMD=sepresultcombine_mTMD;
sepresultcombine.xTMD=sepresultcombine_xTMD;
sepresultcombine.zetaTMD=sepresultcombine_zetaTMD;


a1=evaluate(allresult.prob.Objective,allresult.sol);
a2=evaluate(allresult.prob.Objective,sepresultcombine);

allresultxTMD=xTMDall;
allresult_phi_result=phi_result;
[~,resultall1]=Optim_Damping_for_1_mode(3,sepresultcombine_mTMD,sepresultcombine_zetaTMD,sepresultcombine_fTMD,sepresultcombine_xTMD,[1 2 3]);
[~,resultall2]=Optim_Damping_for_1_mode(3,allresult.sol.mTMD,allresult.sol.zetaTMD,allresult.sol.fTMD,allresult.sol.xTMD,[1 2 3]);

minDamp_aero1=resultall1.minDamp_aero;
minDamp_aero2=resultall2.minDamp_aero;

modes=abs(modes);
sepresultcombine_phi_result=abs(sepresultcombine_phi_result);
allresult_phi_result=abs(allresult_phi_result);

save compare_3tmds_1mode_3modes modes sepresultcombine_xTMD sepresultcombine_phi_result allresultxTMD allresult_phi_result
disp(minDamp_aero1)
disp(minDamp_aero2)
%% 数据分析对比单模态，哪个TMD的安装导致第一阶模态阻尼比骤降

load modeinfo.mat
allresult=importdata("optim_3tmd_123mode_all.mat");

allresultcombine_fTMD=[allresult.sol.fTMD];
allresultcombine_mTMD=[allresult.sol.mTMD];
allresultcombine_xTMD=[allresult.sol.xTMD];
allresultcombine_zetaTMD=[allresult.sol.zetaTMD];

% k1=2;
% allresultcombine_fTMD(k1)=[];
% allresultcombine_mTMD(k1)=[];
% allresultcombine_xTMD(k1)=[];
% allresultcombine_zetaTMD(k1)=[];
[~,resultall]=Optim_Damping_for_1_mode(3,allresultcombine_mTMD,allresultcombine_zetaTMD,allresultcombine_fTMD,allresultcombine_xTMD,[1 2 3]);


sepresult=importdata("optim_1tmd_123mode_separate.mat")
% 模态1
xTMD=sepresult.sol1.xTMD;
[~,index]=sort(abs(nodegap-xTMD));%查找与xTMD最接近的点的排序
xResult=nodegap(index(1:2));%获取最接近的两个点的x坐标
mode2nodes=mode(index(1:2),1:1);%获取两个点坐标的y值
phi_result=interp1(xResult,mode2nodes,xTMD,'linear','extrap');%插值以后任意点的振型

phi_result1=phi_result;
% 模态2
xTMD=sepresult.sol2.xTMD;
[~,index]=sort(abs(nodegap-xTMD));%查找与xTMD最接近的点的排序
xResult=nodegap(index(1:2));%获取最接近的两个点的x坐标
mode2nodes=mode(index(1:2),2);%获取两个点坐标的y值
phi_result=interp1(xResult,mode2nodes,xTMD,'linear','extrap');%插值以后任意点的振型

phi_result2=phi_result;
% 模态3
xTMD=sepresult.sol3.xTMD;
[~,index]=sort(abs(nodegap-xTMD));%查找与xTMD最接近的点的排序
xResult=nodegap(index(1:2));%获取最接近的两个点的x坐标
mode2nodes=mode(index(1:2),3);%获取两个点坐标的y值
phi_result=interp1(xResult,mode2nodes,xTMD,'linear','extrap');%插值以后任意点的振型

phi_result3=phi_result;

sepresultcombine_fTMD=[sepresult.sol1.fTMD sepresult.sol2.fTMD sepresult.sol3.fTMD];
sepresultcombine_mTMD=[sepresult.sol1.mTMD sepresult.sol2.mTMD sepresult.sol3.mTMD];
sepresultcombine_xTMD=[sepresult.sol1.xTMD sepresult.sol2.xTMD sepresult.sol3.xTMD];
sepresultcombine_zetaTMD=[sepresult.sol1.zetaTMD sepresult.sol2.zetaTMD sepresult.sol3.zetaTMD];

% sepresultcombine_fTMD=[sepresult.sol1.fTMD sepresult.sol3.fTMD];
% sepresultcombine_mTMD=[sepresult.sol1.mTMD sepresult.sol3.mTMD];
% sepresultcombine_xTMD=[sepresult.sol1.xTMD sepresult.sol3.xTMD];
% sepresultcombine_zetaTMD=[sepresult.sol1.zetaTMD sepresult.sol3.zetaTMD];


[~,resultall]=Optim_Damping_for_1_mode(3,sepresultcombine_mTMD,sepresultcombine_zetaTMD,sepresultcombine_fTMD,sepresultcombine_xTMD,[1 2 3]);


%% 优化问题(3模态3个TMD)
clc
clear 
close all
if 1
mode_numbers=[1 2 3];%气动力施加的模态，输入n个数，表示分别计算n阶模态
numberofTMD=3;
mass_six_span = 10007779.7;
prob = optimproblem('Description','Optimize the TMDs','ObjectiveSense','maximize'); %优化为最优阻尼比
zetaTMD =optimvar('zetaTMD',numberofTMD,'LowerBound',0.05,'UpperBound',0.25); %阻尼比
fTMD =optimvar('fTMD',numberofTMD,'LowerBound',0.7,'UpperBound',1.5); %频率
mTMD =optimvar('mTMD',numberofTMD,'LowerBound',0.001*mass_six_span,'UpperBound',0.015*mass_six_span); %质量
xTMD = optimvar('xTMD', numberofTMD,'LowerBound', 0, 'UpperBound', 660); %TMD所安装的节点
calmodes_all=3;
% options = optimoptions('particleswarm',Display='iter',PlotFcn='pswplotbestf');
options = optimoptions('particleswarm',Display='iter',PlotFcn='pswplotbestf',UseParallel=true);
% options = optimoptions('ga', 'Display', 'iter', 'PlotFcn', {'gaplotscorediversity', 'gaplotbestf', 'gaplotrange'}, 'UseParallel', true);

Optim_Damping_for_n_foces_n_modes(mode_numbers,1,1,0.01,8.3,50,3)
[minDamping_allmodes,~] = fcn2optimexpr(@Optim_Damping_for_n_foces_n_modes,mode_numbers,numberofTMD,mTMD,zetaTMD,fTMD,xTMD,calmodes_all);


prob.Objective=minDamping_allmodes;
show(prob)

% x0.zetaTMD = 0.1*ones(numberofTMD,1);
% x0.fTMD = 0.82*ones(numberofTMD,1);
% x0.mTMD = 0.015*mass_six_span*ones(numberofTMD,1);
% x0.xTMD = 0.1*ones(numberofTMD,1);

[sol, optval] = solve(prob,'Solver','particleswarm','Options',options);

val = evaluate(prob.Objective,sol);
disp(val)

strmdoes=(num2str(calmodes_all));
strmdoes=strjoin(strsplit(strmdoes),'_');

forcemdoes=(num2str(mode_numbers));
forcemdoes=strjoin(strsplit(forcemdoes),'_');


str="save opt_"+num2str(numberofTMD)+"TMD_"+strmdoes+"modes_forcemodes_"+forcemdoes+".mat";
eval(str)
% save opt_3tmds_3modes
end

%% 数据分析

clc
clear 
close all
load opt_3TMD_3modes_forcemodes_1_2_3.mat

mode_numbers=[1 2 3];
numberofTMD=3;
mTMD=sol.mTMD;
mTMD=ones(3,1)*mTMD(1)*1;
zetaTMD=sol.zetaTMD;
fTMD=sol.fTMD;
xTMD=sol.xTMD;
calmodes_all=3;
[minDamping_allmodes,result]=Optim_Damping_for_n_foces_n_modes(mode_numbers,numberofTMD,mTMD,zetaTMD,fTMD,xTMD,calmodes_all);


load modeinfo.mat

figure
plot(modes(:,1),modes(:,3))
hold on 
plot(modes(:,1),modes(:,4))
plot(modes(:,1),modes(:,5))

[fTMD_sort,sort_seq]=sort(fTMD);
xTMD_sort=xTMD(sort_seq);
mTMD_sort=mTMD(sort_seq);
zetaTMD_sort=zetaTMD(sort_seq);
for k1=1:3
    xTMD_temp=xTMD_sort(k1);
    [~,index]=sort(abs(nodegap-xTMD_temp));%查找与xTMD最接近的点的排序
    xResult=nodegap(index(1:2));%获取最接近的两个点的x坐标
    mode2nodes=mode(index(1:2),k1);%获取两个点坐标的y值
    phi_result(k1,:)=interp1(xResult,mode2nodes,xTMD_temp);%插值以后任意点的振型
end

scatter(xTMD,phi_result)


damping_mode1=result.resultall_modes.result1.result1.Mode.("Damping ratio")
damping_mode2=result.resultall_modes.result2.result1.Mode.("Damping ratio")
damping_mode3=result.resultall_modes.result3.result1.Mode.("Damping ratio")

bardata=[damping_mode1 damping_mode2 damping_mode3]
bardata=bardata'
figure
bar(bardata)


%% 优化问题(3模态3个TMD)-修改为贝叶斯优化
clc
clear 
close all
if 1
mode_numbers=[1 2 3];%气动力施加的模态，输入n个数，表示分别计算n阶模态
numberofTMD=3;
mass_six_span = 10007779.7;
mu=0.015;
mTMD1 = optimizableVariable('mTMD1',[0.001*mass_six_span 0.015*mass_six_span],'Type','real'); %质量
mTMD2 = optimizableVariable('mTMD2',[0.001*mass_six_span 0.015*mass_six_span],'Type','real'); %质量
mTMD3 = optimizableVariable('mTMD3',[0.001*mass_six_span 0.015*mass_six_span],'Type','real'); %质量
zetaTMD1 = optimizableVariable('zetaTMD1',[0.05 0.25],'Type','real'); %阻尼比
zetaTMD2 = optimizableVariable('zetaTMD2',[0.05 0.25],'Type','real'); %阻尼比
zetaTMD3 = optimizableVariable('zetaTMD3',[0.05 0.25],'Type','real'); %阻尼比
fTMD1 = optimizableVariable('fTMD1',[0.7 1.5],'Type','real'); %频率
fTMD2 = optimizableVariable('fTMD2',[0.7 1.5],'Type','real'); %频率
fTMD3 = optimizableVariable('fTMD3',[0.7 1.5],'Type','real'); %频率
xTMD1 = optimizableVariable('xTMD1',[0 660],'Type','real'); %TMD所安装的节点
xTMD2 = optimizableVariable('xTMD2',[0 660],'Type','real'); %TMD所安装的节点
xTMD3 = optimizableVariable('xTMD3',[0 660],'Type','real'); %TMD所安装的节点
calmodes_all=3;
% mTMD=[mTMD1 mTMD2 mTMD3];
% zetaTMD=[zetaTMD1 zetaTMD2 zetaTMD3];
% fTMD=[fTMD1 fTMD2 fTMD3];
% xTMD=[xTMD1 xTMD2 xTMD3];

fun=@(x)Optim_Damping_for_n_foces_n_modes_bayesopt(mode_numbers,numberofTMD,x.mTMD1,x.mTMD2,x.mTMD3,x.zetaTMD1,x.zetaTMD2,x.zetaTMD3,x.fTMD1,x.fTMD2,x.fTMD3,x.xTMD1,x.xTMD2,x.xTMD3,calmodes_all);
% fun2=@(X)xconstraint(X,mu,mass_six_span);
results = bayesopt(fun,[mTMD1 mTMD2 mTMD3 zetaTMD1 zetaTMD2 zetaTMD3 fTMD1 fTMD2 fTMD3 xTMD1 xTMD2 xTMD3],'XConstraintFcn',@xconstraint,'AcquisitionFunctionName','expected-improvement-plus','MaxObjectiveEvaluations',100,'UseParallel',true);





% prob = optimproblem('Description','Optimize the TMDs','ObjectiveSense','maximize'); %优化为最优阻尼比
% zetaTMD =optimvar('zetaTMD',numberofTMD,'LowerBound',0.05,'UpperBound',0.25); %阻尼比
% fTMD =optimvar('fTMD',numberofTMD,'LowerBound',0.7,'UpperBound',1.5); %频率
% mTMD =optimvar('mTMD',numberofTMD,'LowerBound',0.001*mass_six_span,'UpperBound',0.015*mass_six_span); %质量
% xTMD = optimvar('xTMD', numberofTMD,'LowerBound', 0, 'UpperBound', 660); %TMD所安装的节点
% calmodes_all=3;
% % options = optimoptions('particleswarm',Display='iter',PlotFcn='pswplotbestf');
% options = optimoptions('particleswarm',Display='iter',PlotFcn='pswplotbestf',UseParallel=true);
% % options = optimoptions('ga', 'Display', 'iter', 'PlotFcn', {'gaplotscorediversity', 'gaplotbestf', 'gaplotrange'}, 'UseParallel', true);

% Optim_Damping_for_n_foces_n_modes(mode_numbers,1,1,0.01,8.3,50,3)
% [minDamping_allmodes,~] = fcn2optimexpr(@Optim_Damping_for_n_foces_n_modes,mode_numbers,numberofTMD,mTMD,zetaTMD,fTMD,xTMD,calmodes_all);


% prob.Objective=minDamping_allmodes;
% show(prob)

% % x0.zetaTMD = 0.1*ones(numberofTMD,1);
% % x0.fTMD = 0.82*ones(numberofTMD,1);
% % x0.mTMD = 0.015*mass_six_span*ones(numberofTMD,1);
% % x0.xTMD = 0.1*ones(numberofTMD,1);

% [sol, optval] = solve(prob,'Solver','particleswarm','Options',options);

% val = evaluate(prob.Objective,sol);
% disp(val)

% strmdoes=(num2str(calmodes_all));
% strmdoes=strjoin(strsplit(strmdoes),'_');

% forcemdoes=(num2str(mode_numbers));
% forcemdoes=strjoin(strsplit(forcemdoes),'_');


% str="save opt_"+num2str(numberofTMD)+"TMD_"+strmdoes+"modes_forcemodes_"+forcemdoes+".mat";
% eval(str)
% save opt_3tmds_3modes
end



%% 所用到的函数
function tf = xconstraint(x)
    totalmass=x.mTMD1+x.mTMD2+x.mTMD3;
    mu=0.015;
    mass_six_span = 10007779.7;
    tf = totalmass<=mu*mass_six_span & totalmass>=mu*0.95*mass_six_span;
end