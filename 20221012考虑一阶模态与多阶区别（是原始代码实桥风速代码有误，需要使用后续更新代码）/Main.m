%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%Date: 2022-10-13 11:11:22
%LastEditors: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%LastEditTime: 2022-10-13 19:57:43
%FilePath: \NonlinearScanlan\20221012考虑一阶模态与多阶区别\Main.m
%Description: 主文件（调用对比函数进行批量分析）
%
%Copyright (c) 2022 by xushengyichn 54436848+xushengyichn@users.noreply.github.com, All Rights Reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%% 计算两阶模态 安装于1036节点
clc;clear;close all;
addpath("..\函数\")
calmodes=2;
mus=0.01/100:0.01/100:5/100;
Frequency=zeros(calmodes+1,length(mus));
Damping_ratio=zeros(calmodes+1,length(mus));
dis=zeros(1,length(mus));
% for k1 = 1:1

numIterations = length(mus);
res = zeros(numIterations, 1);


ppm = ParforProgressbar(numIterations,'showWorkerProgress',true,'progressBarUpdatePeriod',3,'title','my fancy title'); 
pauseTime = 60/numIterations;

parfor k1 = 1:numIterations
    mu=mus(k1);
    zetaTMD=0.06;
    fTMD=0.8339;
    locationTMD=1036;
    t_length=500;
    flag = 0;%判断是否收敛
    result=[];
    disnode_sum=[];
    dis1=[];
    dis2=[];
    while flag==0
    result=Compare_n_modes(mu,zetaTMD,fTMD,locationTMD,calmodes,t_length);
    disnode=result.dis_nodeTMD;
    disnode_1=disnode(1,:);
    disnode_2=disnode(2,:);
    disnode_sum=disnode_1+disnode_2;
    dis1=max(abs(disnode_sum(round(end-200,0):round(end-100,0))));
    dis2=max(abs(disnode_sum(round(end-100,0):end)));
    decided_value=abs(dis1/dis2);
    if and( decided_value>0.95, decided_value<1.05)
        flag=1;
         disp(num2str(mu)+"已收敛");
    else
        t_length=t_length*2;
        disp(num2str(mu)+"未收敛，计算时间调整为："+num2str(t_length));
    end
    end
    Frequency(:,k1)=result.output.Mode.Frequency;
    Damping_ratio(:,k1)=result.output.Mode.("Damping ratio");
    dis(k1)=dis2;
        pause(pauseTime);
    % increment counter to track progress
    ppm.increment();
end

collectresult=[mus' dis' Frequency' Damping_ratio'];
resulttable=array2table(collectresult,"VariableNames",{'mass ratio','displacement','fre1','fre2','fre3','zeta1','zeta2','zeta3'});
save TMD1036_2mode_chageMassRatio resulttable

%% 只计算一阶模态 安装于1036节点

clc;clear;close all;
addpath("..\函数\")
calmodes=1;
mus=0.01/100:0.01/100:5/100;
Frequency=zeros(calmodes+1,length(mus));
Damping_ratio=zeros(calmodes+1,length(mus));
dis=zeros(1,length(mus));
% for k1 = 1:1

numIterations = length(mus);
res = zeros(numIterations, 1);


ppm = ParforProgressbar(numIterations,'showWorkerProgress',true,'progressBarUpdatePeriod',3,'title','my fancy title'); 
pauseTime = 60/numIterations;
for k1 = 100
% parfor k1 = 1:numIterations
    mus=0.01/100:0.01/100:5/100;
    mu=mus(k1);
    zetaTMD=0.06;
    fTMD=0.8339;
    locationTMD=1036;
    t_length=100;
    flag = 0;%判断是否收敛
    result=[];
    disnode_sum=[];
    dis1=[];
    dis2=[];
    while flag==0
    result=Compare_n_modes(mu,zetaTMD,fTMD,locationTMD,calmodes,t_length);
    disnode=result.dis_nodeTMD;
    disnode_1=disnode(1,:);
    disnode_sum=disnode_1;
    dis1=max(abs(disnode_sum(round(end-200,0):round(end-100,0))));
    dis2=max(abs(disnode_sum(round(end-100,0):end)));
    decided_value=abs(dis1/dis2);
        if and( decided_value>0.95, decided_value<1.05)
            flag=1;
             disp(num2str(mu)+"已收敛");
        else
            t_length=t_length*2;
            disp(num2str(mu)+"未收敛，计算时间调整为："+num2str(t_length));
        end
    end
    Frequency(:,k1)=result.output.Mode.Frequency;
    Damping_ratio(:,k1)=result.output.Mode.("Damping ratio");
    dis(k1)=dis2;
        pause(pauseTime);
    % increment counter to track progress
    ppm.increment();
end

collectresult=[mus' dis' Frequency' Damping_ratio'];
resulttable=array2table(collectresult,"VariableNames",{'mass ratio','displacement','fre1','fre2','zeta1','zeta2'});
save TMD1036_1mode_chageMassRatio resulttable



%% 计算两阶模态 安装于1162节点
clc;clear;close all;
addpath("..\函数\")
calmodes=2;
mus=0.01/100:0.01/100:5/100;
Frequency=zeros(calmodes+1,length(mus));
Damping_ratio=zeros(calmodes+1,length(mus));
dis=zeros(1,length(mus));
% for k1 = 1:1

numIterations = length(mus);
res = zeros(numIterations, 1);


ppm = ParforProgressbar(numIterations,'showWorkerProgress',true,'progressBarUpdatePeriod',3,'title','my fancy title'); 
pauseTime = 60/numIterations;

parfor k1 = 1:numIterations
    mu=mus(k1);
    zetaTMD=0.06;
    fTMD=0.8339;
    locationTMD=1162;
    t_length=100;
    flag = 0;%判断是否收敛
    result=[];
    disnode_sum=[];
    dis1=[];
    dis2=[];
    while flag==0
    result=Compare_n_modes(mu,zetaTMD,fTMD,locationTMD,calmodes,t_length);
    disnode=result.dis_nodeTMD;
    disnode_1=disnode(1,:);
    disnode_2=disnode(2,:);
    disnode_sum=disnode_1+disnode_2;
    dis1=max(abs(disnode_sum(round(end-200,0):round(end-100,0))));
    dis2=max(abs(disnode_sum(round(end-100,0):end)));
    decided_value=abs(dis1/dis2);
        if and( decided_value>0.95, decided_value<1.05)
            flag=1;
             disp(num2str(mu)+"已收敛");
        else
            t_length=t_length*2;
            disp(num2str(mu)+"未收敛，计算时间调整为："+num2str(t_length));
        end
    end
    Frequency(:,k1)=result.output.Mode.Frequency;
    Damping_ratio(:,k1)=result.output.Mode.("Damping ratio");
    dis(k1)=dis2;
        pause(pauseTime);
    % increment counter to track progress
    ppm.increment();
end

collectresult=[mus' dis' Frequency' Damping_ratio'];
resulttable=array2table(collectresult,"VariableNames",{'mass ratio','displacement','fre1','fre2','fre3','zeta1','zeta2','zeta3'});
save TMD1162_2mode_chageMassRatio resulttable

%% 只计算一阶模态 安装于1162节点

clc;clear;close all;
addpath("..\函数\")
calmodes=1;
mus=0.01/100:0.01/100:5/100;
Frequency=zeros(calmodes+1,length(mus));
Damping_ratio=zeros(calmodes+1,length(mus));
dis=zeros(1,length(mus));
% for k1 = 1:1

numIterations = length(mus);
res = zeros(numIterations, 1);


ppm = ParforProgressbar(numIterations,'showWorkerProgress',true,'progressBarUpdatePeriod',3,'title','my fancy title'); 
pauseTime = 60/numIterations;

parfor k1 = 1:numIterations
    mu=mus(k1);
    zetaTMD=0.06;
    fTMD=0.8339;
    locationTMD=1162;
    t_length=100;
    flag = 0;%判断是否收敛
    result=[];
    disnode_sum=[];
    dis1=[];
    dis2=[];
    while flag==0
    result=Compare_n_modes(mu,zetaTMD,fTMD,locationTMD,calmodes,t_length);
    disnode=result.dis_nodeTMD;
    disnode_1=disnode(1,:);
    disnode_sum=disnode_1;
    dis1=max(abs(disnode_sum(round(end-200,0):round(end-100,0))));
    dis2=max(abs(disnode_sum(round(end-100,0):end)));
    decided_value=abs(dis1/dis2);
        if and( decided_value>0.95, decided_value<1.05)
            flag=1;
             disp(num2str(mu)+"已收敛");
        else
            t_length=t_length*2;
            disp(num2str(mu)+"未收敛，计算时间调整为："+num2str(t_length));
        end
    end
    Frequency(:,k1)=result.output.Mode.Frequency;
    Damping_ratio(:,k1)=result.output.Mode.("Damping ratio");
    dis(k1)=dis2;
        pause(pauseTime);
    % increment counter to track progress
    ppm.increment();
end

collectresult=[mus' dis' Frequency' Damping_ratio'];
resulttable=array2table(collectresult,"VariableNames",{'mass ratio','displacement','fre1','fre2','zeta1','zeta2'});
save TMD1162_1mode_chageMassRatio resulttable


%% 计算两阶模态 安装于1179节点
clc;clear;close all;
addpath("..\函数\")
calmodes=2;
mus=0.01/100:0.01/100:5/100;
Frequency=zeros(calmodes+1,length(mus));
Damping_ratio=zeros(calmodes+1,length(mus));
dis=zeros(1,length(mus));
% for k1 = 1:1

numIterations = length(mus);
res = zeros(numIterations, 1);


ppm = ParforProgressbar(numIterations,'showWorkerProgress',true,'progressBarUpdatePeriod',3,'title','my fancy title'); 
pauseTime = 60/numIterations;

parfor k1 = 1:numIterations
    mu=mus(k1);
    zetaTMD=0.06;
    fTMD=0.8339;
    locationTMD=1179;
    t_length=100;
    flag = 0;%判断是否收敛
    result=[];
    disnode_sum=[];
    dis1=[];
    dis2=[];
    while flag==0
    result=Compare_n_modes(mu,zetaTMD,fTMD,locationTMD,calmodes,t_length);
    disnode=result.dis_nodeTMD;
    disnode_1=disnode(1,:);
    disnode_2=disnode(2,:);
    disnode_sum=disnode_1+disnode_2;
    dis1=max(abs(disnode_sum(round(end-200,0):round(end-100,0))));
    dis2=max(abs(disnode_sum(round(end-100,0):end)));
    decided_value=abs(dis1/dis2);
        if and( decided_value>0.95, decided_value<1.05)
            flag=1;
             disp(num2str(mu)+"已收敛");
        else
            t_length=t_length*2;
            disp(num2str(mu)+"未收敛，计算时间调整为："+num2str(t_length));
        end
    end
    Frequency(:,k1)=result.output.Mode.Frequency;
    Damping_ratio(:,k1)=result.output.Mode.("Damping ratio");
    dis(k1)=dis2;
        pause(pauseTime);
    % increment counter to track progress
    ppm.increment();
end

collectresult=[mus' dis' Frequency' Damping_ratio'];
resulttable=array2table(collectresult,"VariableNames",{'mass ratio','displacement','fre1','fre2','fre3','zeta1','zeta2','zeta3'});
save TMD1179_2mode_chageMassRatio resulttable


%% 只计算一阶模态 安装于1179节点

clc;clear;close all;
addpath("..\函数\")
calmodes=1;
mus=0.01/100:0.01/100:5/100;
Frequency=zeros(calmodes+1,length(mus));
Damping_ratio=zeros(calmodes+1,length(mus));
dis=zeros(1,length(mus));
% for k1 = 1:1

numIterations = length(mus);
res = zeros(numIterations, 1);


ppm = ParforProgressbar(numIterations,'showWorkerProgress',true,'progressBarUpdatePeriod',3,'title','my fancy title'); 
pauseTime = 60/numIterations;

parfor k1 = 1:numIterations
    mu=mus(k1);
    zetaTMD=0.06;
    fTMD=0.8339;
    locationTMD=1179;
    t_length=100;
    flag = 0;%判断是否收敛
    result=[];
    disnode_sum=[];
    dis1=[];
    dis2=[];
    while flag==0
    result=Compare_n_modes(mu,zetaTMD,fTMD,locationTMD,calmodes,t_length);
    disnode=result.dis_nodeTMD;
    disnode_1=disnode(1,:);
    disnode_sum=disnode_1;
    dis1=max(abs(disnode_sum(round(end-200,0):round(end-100,0))));
    dis2=max(abs(disnode_sum(round(end-100,0):end)));
    decided_value=abs(dis1/dis2);
        if and( decided_value>0.95, decided_value<1.05)
            flag=1;
             disp(num2str(mu)+"已收敛");
        else
            t_length=t_length*2;
            disp(num2str(mu)+"未收敛，计算时间调整为："+num2str(t_length));
        end
    end
    Frequency(:,k1)=result.output.Mode.Frequency;
    Damping_ratio(:,k1)=result.output.Mode.("Damping ratio");
    dis(k1)=dis2;
        pause(pauseTime);
    % increment counter to track progress
    ppm.increment();
end

collectresult=[mus' dis' Frequency' Damping_ratio'];
resulttable=array2table(collectresult,"VariableNames",{'mass ratio','displacement','fre1','fre2','zeta1','zeta2'});
save TMD1179_1mode_chageMassRatio resulttable


%% 计算两阶模态 安装于1036节点 TMD参数会变动
clc;clear;close all;
addpath("..\函数\")
calmodes=2;
mus=0.01/100:0.01/100:5/100;
Frequency=zeros(calmodes+1,length(mus));
Damping_ratio=zeros(calmodes+1,length(mus));
dis=zeros(1,length(mus));
% for k1 = 1:1

numIterations = length(mus);
res = zeros(numIterations, 1);

modedata=importdata("modedata.mat");
ppm = ParforProgressbar(numIterations,'showWorkerProgress',true,'progressBarUpdatePeriod',3,'title','my fancy title'); 
pauseTime = 60/numIterations;

% parfor k1 = 100
parfor k1 = 1:numIterations
    mu=mus(k1);
    locationTMD=1036;
    phiTMD=modedata(find(locationTMD==modedata(:,1)),3);
    muTMD_real=10007779.7*mu*phiTMD^2;
    zetaTMD=sqrt(3*muTMD_real/8/(1+muTMD_real));
    fTMD=0.8339*(1+muTMD_real);
    t_length=500;
    flag = 0;%判断是否收敛
    result=[];
    disnode_sum=[];
    dis1=[];
    dis2=[];
    while flag==0
    result=Compare_n_modes(mu,zetaTMD,fTMD,locationTMD,calmodes,t_length);
    disnode=result.dis_nodeTMD;
    disnode_1=disnode(1,:);
    disnode_2=disnode(2,:);
    disnode_sum=disnode_1+disnode_2;
    dis1=max(abs(disnode_sum(round(end-200,0):round(end-100,0))));
    dis2=max(abs(disnode_sum(round(end-100,0):end)));
    decided_value=abs(dis1/dis2);
    if and( decided_value>0.95, decided_value<1.05)
        flag=1;
         disp(num2str(mu)+"已收敛");
    else
        t_length=t_length*2;
        disp(num2str(mu)+"未收敛，计算时间调整为："+num2str(t_length));
    end
    end
    Frequency(:,k1)=result.output.Mode.Frequency;
    Damping_ratio(:,k1)=result.output.Mode.("Damping ratio");
    dis(k1)=dis2;
        pause(pauseTime);
    % increment counter to track progress
    ppm.increment();
end

collectresult=[mus' dis' Frequency' Damping_ratio'];
resulttable=array2table(collectresult,"VariableNames",{'mass ratio','displacement','fre1','fre2','fre3','zeta1','zeta2','zeta3'});
save TMD1036_2mode_chageMassRatio resulttable

%% 只计算一阶模态 安装于1036节点 TMD参数会变动

clc;clear;close all;
addpath("..\函数\")
calmodes=1;
mus=0.01/100:0.01/100:5/100;
Frequency=zeros(calmodes+1,length(mus));
Damping_ratio=zeros(calmodes+1,length(mus));
dis=zeros(1,length(mus));
% for k1 = 1:1

numIterations = length(mus);
res = zeros(numIterations, 1);

modedata=importdata("modedata.mat");
ppm = ParforProgressbar(numIterations,'showWorkerProgress',true,'progressBarUpdatePeriod',3,'title','my fancy title'); 
pauseTime = 60/numIterations;

for k1 = 100
% parfor k1 = 1:numIterations
     mu=mus(k1);
     mu=0.01
    locationTMD=1036;
    phiTMD=modedata(find(locationTMD==modedata(:,1)),3);
    muTMD_real=10007779.7*mu*phiTMD^2;
    zetaTMD=sqrt(3*muTMD_real/8/(1+muTMD_real));
    fTMD=0.8339*(1+muTMD_real);
    t_length=500;
    flag = 0;%判断是否收敛
    result=[];
    disnode_sum=[];
    dis1=[];
    dis2=[];
    while flag==0
    result=Compare_n_modes(mu,zetaTMD,fTMD,locationTMD,calmodes,t_length);
    disnode=result.dis_nodeTMD;
    disnode_1=disnode(1,:);
    disnode_sum=disnode_1;
    dis1=max(abs(disnode_sum(round(end-200,0):round(end-100,0))));
    dis2=max(abs(disnode_sum(round(end-100,0):end)));
    decided_value=abs(dis1/dis2);
        if and( decided_value>0.95, decided_value<1.05)
            flag=1;
             disp(num2str(mu)+"已收敛");
        else
            t_length=t_length*2;
            disp(num2str(mu)+"未收敛，计算时间调整为："+num2str(t_length));
        end
    end
    Frequency(:,k1)=result.output.Mode.Frequency;
    Damping_ratio(:,k1)=result.output.Mode.("Damping ratio");
    dis(k1)=dis2;
        pause(pauseTime);
    % increment counter to track progress
    ppm.increment();
end

collectresult=[mus' dis' Frequency' Damping_ratio'];
resulttable=array2table(collectresult,"VariableNames",{'mass ratio','displacement','fre1','fre2','zeta1','zeta2'});
save TMD1036_1mode_chageMassRatio resulttable



%% 计算两阶模态 安装于1162节点
clc;clear;close all;
addpath("..\函数\")
calmodes=2;
mus=0.01/100:0.01/100:5/100;
Frequency=zeros(calmodes+1,length(mus));
Damping_ratio=zeros(calmodes+1,length(mus));
dis=zeros(1,length(mus));
% for k1 = 1:1

numIterations = length(mus);
res = zeros(numIterations, 1);


ppm = ParforProgressbar(numIterations,'showWorkerProgress',true,'progressBarUpdatePeriod',3,'title','my fancy title'); 
pauseTime = 60/numIterations;

parfor k1 = 1:numIterations
    mu=mus(k1);
    zetaTMD=0.06;
    fTMD=0.8339;
    locationTMD=1162;
    t_length=100;
    flag = 0;%判断是否收敛
    result=[];
    disnode_sum=[];
    dis1=[];
    dis2=[];
    while flag==0
    result=Compare_n_modes(mu,zetaTMD,fTMD,locationTMD,calmodes,t_length);
    disnode=result.dis_nodeTMD;
    disnode_1=disnode(1,:);
    disnode_2=disnode(2,:);
    disnode_sum=disnode_1+disnode_2;
    dis1=max(abs(disnode_sum(round(end-200,0):round(end-100,0))));
    dis2=max(abs(disnode_sum(round(end-100,0):end)));
    decided_value=abs(dis1/dis2);
        if and( decided_value>0.95, decided_value<1.05)
            flag=1;
             disp(num2str(mu)+"已收敛");
        else
            t_length=t_length*2;
            disp(num2str(mu)+"未收敛，计算时间调整为："+num2str(t_length));
        end
    end
    Frequency(:,k1)=result.output.Mode.Frequency;
    Damping_ratio(:,k1)=result.output.Mode.("Damping ratio");
    dis(k1)=dis2;
        pause(pauseTime);
    % increment counter to track progress
    ppm.increment();
end

collectresult=[mus' dis' Frequency' Damping_ratio'];
resulttable=array2table(collectresult,"VariableNames",{'mass ratio','displacement','fre1','fre2','fre3','zeta1','zeta2','zeta3'});
save TMD1162_2mode_chageMassRatio resulttable

%% 只计算一阶模态 安装于1162节点

clc;clear;close all;
addpath("..\函数\")
calmodes=1;
mus=0.01/100:0.01/100:5/100;
Frequency=zeros(calmodes+1,length(mus));
Damping_ratio=zeros(calmodes+1,length(mus));
dis=zeros(1,length(mus));
% for k1 = 1:1

numIterations = length(mus);
res = zeros(numIterations, 1);


ppm = ParforProgressbar(numIterations,'showWorkerProgress',true,'progressBarUpdatePeriod',3,'title','my fancy title'); 
pauseTime = 60/numIterations;

parfor k1 = 1:numIterations
    mu=mus(k1);
    zetaTMD=0.06;
    fTMD=0.8339;
    locationTMD=1162;
    t_length=100;
    flag = 0;%判断是否收敛
    result=[];
    disnode_sum=[];
    dis1=[];
    dis2=[];
    while flag==0
    result=Compare_n_modes(mu,zetaTMD,fTMD,locationTMD,calmodes,t_length);
    disnode=result.dis_nodeTMD;
    disnode_1=disnode(1,:);
    disnode_sum=disnode_1;
    dis1=max(abs(disnode_sum(round(end-200,0):round(end-100,0))));
    dis2=max(abs(disnode_sum(round(end-100,0):end)));
    decided_value=abs(dis1/dis2);
        if and( decided_value>0.95, decided_value<1.05)
            flag=1;
             disp(num2str(mu)+"已收敛");
        else
            t_length=t_length*2;
            disp(num2str(mu)+"未收敛，计算时间调整为："+num2str(t_length));
        end
    end
    Frequency(:,k1)=result.output.Mode.Frequency;
    Damping_ratio(:,k1)=result.output.Mode.("Damping ratio");
    dis(k1)=dis2;
        pause(pauseTime);
    % increment counter to track progress
    ppm.increment();
end

collectresult=[mus' dis' Frequency' Damping_ratio'];
resulttable=array2table(collectresult,"VariableNames",{'mass ratio','displacement','fre1','fre2','zeta1','zeta2'});
save TMD1162_1mode_chageMassRatio resulttable


%% 计算两阶模态 安装于1179节点
clc;clear;close all;
addpath("..\函数\")
calmodes=2;
mus=0.01/100:0.01/100:5/100;
Frequency=zeros(calmodes+1,length(mus));
Damping_ratio=zeros(calmodes+1,length(mus));
dis=zeros(1,length(mus));
% for k1 = 1:1

numIterations = length(mus);
res = zeros(numIterations, 1);


ppm = ParforProgressbar(numIterations,'showWorkerProgress',true,'progressBarUpdatePeriod',3,'title','my fancy title'); 
pauseTime = 60/numIterations;

parfor k1 = 1:numIterations
    mu=mus(k1);
    zetaTMD=0.06;
    fTMD=0.8339;
    locationTMD=1179;
    t_length=100;
    flag = 0;%判断是否收敛
    result=[];
    disnode_sum=[];
    dis1=[];
    dis2=[];
    while flag==0
    result=Compare_n_modes(mu,zetaTMD,fTMD,locationTMD,calmodes,t_length);
    disnode=result.dis_nodeTMD;
    disnode_1=disnode(1,:);
    disnode_2=disnode(2,:);
    disnode_sum=disnode_1+disnode_2;
    dis1=max(abs(disnode_sum(round(end-200,0):round(end-100,0))));
    dis2=max(abs(disnode_sum(round(end-100,0):end)));
    decided_value=abs(dis1/dis2);
        if and( decided_value>0.95, decided_value<1.05)
            flag=1;
             disp(num2str(mu)+"已收敛");
        else
            t_length=t_length*2;
            disp(num2str(mu)+"未收敛，计算时间调整为："+num2str(t_length));
        end
    end
    Frequency(:,k1)=result.output.Mode.Frequency;
    Damping_ratio(:,k1)=result.output.Mode.("Damping ratio");
    dis(k1)=dis2;
        pause(pauseTime);
    % increment counter to track progress
    ppm.increment();
end

collectresult=[mus' dis' Frequency' Damping_ratio'];
resulttable=array2table(collectresult,"VariableNames",{'mass ratio','displacement','fre1','fre2','fre3','zeta1','zeta2','zeta3'});
save TMD1179_2mode_chageMassRatio resulttable


%% 只计算一阶模态 安装于1179节点

clc;clear;close all;
addpath("..\函数\")
calmodes=1;
mus=0.01/100:0.01/100:5/100;
Frequency=zeros(calmodes+1,length(mus));
Damping_ratio=zeros(calmodes+1,length(mus));
dis=zeros(1,length(mus));
% for k1 = 1:1

numIterations = length(mus);
res = zeros(numIterations, 1);


ppm = ParforProgressbar(numIterations,'showWorkerProgress',true,'progressBarUpdatePeriod',3,'title','my fancy title'); 
pauseTime = 60/numIterations;

parfor k1 = 1:numIterations
    mu=mus(k1);
    zetaTMD=0.06;
    fTMD=0.8339;
    locationTMD=1179;
    t_length=100;
    flag = 0;%判断是否收敛
    result=[];
    disnode_sum=[];
    dis1=[];
    dis2=[];
    while flag==0
    result=Compare_n_modes(mu,zetaTMD,fTMD,locationTMD,calmodes,t_length);
    disnode=result.dis_nodeTMD;
    disnode_1=disnode(1,:);
    disnode_sum=disnode_1;
    dis1=max(abs(disnode_sum(round(end-200,0):round(end-100,0))));
    dis2=max(abs(disnode_sum(round(end-100,0):end)));
    decided_value=abs(dis1/dis2);
        if and( decided_value>0.95, decided_value<1.05)
            flag=1;
             disp(num2str(mu)+"已收敛");
        else
            t_length=t_length*2;
            disp(num2str(mu)+"未收敛，计算时间调整为："+num2str(t_length));
        end
    end
    Frequency(:,k1)=result.output.Mode.Frequency;
    Damping_ratio(:,k1)=result.output.Mode.("Damping ratio");
    dis(k1)=dis2;
        pause(pauseTime);
    % increment counter to track progress
    ppm.increment();
end

collectresult=[mus' dis' Frequency' Damping_ratio'];
resulttable=array2table(collectresult,"VariableNames",{'mass ratio','displacement','fre1','fre2','zeta1','zeta2'});
save TMD1179_1mode_chageMassRatio resulttable