%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%Date: 2022-09-12 10:21:29
%LastEditors: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%LastEditTime: 2022-09-12 10:46:00
%FilePath: \NonlinearScanlan\Complex_Eigenvalue_Analysis.m
%Description: 给定MCK矩阵求特征值，返回考虑阻尼矩阵情况下，每一阶模态的频率和阻尼比
%
%Copyright (c) 2022 by xushengyichn 54436848+xushengyichn@users.noreply.github.com, All Rights Reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Mode = Complex_Eigenvalue_Analysis(M,C,K)

%COMPLEX_EIGENVALUE_ANALYSIS 给定MCK矩阵求特征值
%   M,C,K为矩阵
%   eigenvector为特征向量
%   eigenvalue为特征值


%    %调试代码使用该行
% clc;clear all;close all;
%     load MCK.mat
%     M=MM;
%     C=CC;
%     K=KK;

    matrixsize=size(M);
    A=[zeros(matrixsize) -M;-M -C];%矩阵A
    B=[-M zeros(matrixsize);zeros(matrixsize) K];%矩阵B
    A_1=inv(A);
    % 求解特征值
    A_1B=A_1*B;
    % A_1B2=[-inv(M)*C -inv(M)*K;eye(matrixsize) zeros(matrixsize)]; %与上式相同
    [V, D]=eig(A_1B);
    [sn,ind] = sort(diag(D));
    Ds = diag(D(ind,ind)); %特征值
    Vs = V(:,ind); %特征向量

    % 特征向量归一化
    temp=conj(Vs)'*A*Vs;
    temp2=diag(sqrt(temp));
    for k1 = 1:matrixsize*2
        Vs(:,k1)=Vs(:,k1)./temp2(k1);
    end

    omega=unique(abs(imag(Ds)),'stable');
    Fre=omega/2/pi;
    xi_omega=unique(real(Ds),'stable');
    xi=xi_omega./omega*-1;
    
    modenum=(1:1:matrixsize)';
    mode=[modenum Fre xi];
    TableName={'Mode','Frequency','Damping ratio'};
    Mode=array2table(mode,'VariableNames',TableName);
    disp(Mode)
end

    % 以下为计算自由振动计算初始状态以及自由振动响应的代码，此处不需要
    % x0=[0;0;1;2];%初始状态
    % disp("initial condition")
    % disp(x0)
    % disp("Matrix A")
    % disp(A)

    % d=conj(Vs)'*A*x0;
    % disp("d")
    % disp(d)



    % syms t; % be careful performance issue.

    % sn
    % for k1=1:4
    %     alpha(k1)=real(sn(k1));
    %     beta(k1)=imag(sn(k1));
    % end

    % x1=d(1)*exp(alpha(1)*t)*(cos(beta(1)*t)+sin(beta(1)*t)*1i)*Vs(:,1);
    % x2=d(2)*exp(alpha(2)*t)*(cos(beta(2)*t)+sin(beta(2)*t)*1i)*Vs(:,2);
    % x=x1+x2;
    % x=simple(x);
    % x=vpa(x,4) ;
    % disp("----")
    % disp(x(3))

    % % u1d0t=0;
    % % u2dot=0;
    % % u1=0;
    % % u2=0;
    % % for k1 =1:4
    % %     sum=d(k1)*exp(alpha(k1)*t)*(cos(beta(k1)*t)+sin(beta(k1)*t)*1i)+sum;
    % % end
    % % sum=simplify(sum);
    % % sum=vpa(sum,4) 
    % % disp(sum);


