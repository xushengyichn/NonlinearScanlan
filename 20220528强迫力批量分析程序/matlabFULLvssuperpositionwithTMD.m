% 本代码使用ansys导出矩阵计算振动响应，分别使用FUll矩阵和模态叠加法，计算结果不一致，需要采用复模态的方法，因为阻尼矩阵不是正交的，代码未修改。
clc
clear
close all

hb_to_mm ( 'KMatrixTMD.matrix', 'KTMD.txt' );
hb_to_mm ( 'MMatrixTMD.matrix', 'MTMD.txt' );
hb_to_mm ( 'CMatrixTMD.matrix', 'CTMD.txt' );



Kdata = importdata('KTMD.txt').data;
Kmatrix = zeros(Kdata(1,1),Kdata(1,2));
for i = 2:size(Kdata,1)
    Kmatrix(Kdata(i,1),Kdata(i,2)) = Kdata(i,3);
end

Mdata = importdata('MTMD.txt').data;
Mmatrix = zeros(Mdata(1,1),Mdata(1,2));

for i = 2:size(Mdata,1)
    Mmatrix(Mdata(i,1),Mdata(i,2)) = Mdata(i,3);
end

Cdata = importdata('CTMD.txt').data;
Cmatrix = zeros(Cdata(1,1),Cdata(1,2));

for i = 2:size(Cdata,1)
    Cmatrix(Cdata(i,1),Cdata(i,2)) = Cdata(i,3);
end

% 还原对角线以上元素，使之为对称阵, ANSYS只给出下三角矩阵
K =diag(diag(Kmatrix)/2)+Kmatrix-diag(diag(Kmatrix));
K = K+K';
M =diag(diag(Mmatrix)/2)+Mmatrix-diag(diag(Mmatrix));
M = M+M';
C =diag(diag(Cmatrix)/2)+Cmatrix-diag(diag(Cmatrix));
C = C+C';


fangdaxishu=10;

PI=3.14159265359;
FRE=0.069704171453635;	

THETA=2*PI*FRE;

dt=0.01;		
% !计算时间（秒）
T=300;
NNT=T/dt;
t=0:dt:T;
P1=fangdaxishu*sin(THETA*t);
P2=fangdaxishu*sin(2*THETA*t)*50;
degrees=size(M,1);
P_eachpoint=zeros(degrees,length(t));
P_eachpoint(1,:)=P1;
P_eachpoint(141,:)=P2;
P=P_eachpoint;
% P=zeros(50,length(t));
% for t1=1:matrixsize
%     if t1<=nModes
%         P(t1,:)=Peq(t1,mode_vec3,KMmapping,P_eachpoint,points,t);
%     else
%         P(t1,:)=zeros(1,size(t,2));
%     end
% end
KMmapping = importmappingmatrix('KMatrixTMD.mapping');
u0=zeros(degrees,1);
udot0=zeros(degrees,1);
nbeta=0.25;
ngam=0.5;

[u udot u2dot] = NewmarkInt(t,M,C,K,P,ngam,nbeta,u0,udot0);
clear u0 udot0 P_eachpoint P1

% 特征值分析，即计算频率Freq和振型Phi，50代表求50阶，SM表示从较小的特征值开始求解
calmodes=50;
matrixsize=calmodes;
[eig_vec,eig_val]=eigs(K,M,calmodes,'SM');
[nfdof,nfdof]=size(eig_vec);
for j=1:nfdof
    mnorm=sqrt(eig_vec(:,j)'*M*eig_vec(:,j));
    eig_vec(:,j)=eig_vec(:,j)/mnorm;%振型质量归一化
end
[omeg,w_order] =sort(sqrt(diag(eig_val)));
mode_vec=eig_vec(:,w_order);
Freq = omeg/(2*pi);


%% Create the matrix by superposition method
MM=zeros(matrixsize,matrixsize);
for t1=1:matrixsize
        MM(t1,t1)=P_eq(t1,mode_vec,M);
end
clear t1 


CC=zeros(matrixsize,matrixsize);
for t1=1:matrixsize
        CC(t1,t1)=P_eq(t1,mode_vec,C);
end
clear t1 

KK=zeros(matrixsize,matrixsize);
for t1=1:matrixsize
        KK(t1,t1)=P_eq(t1,mode_vec,K);
end
clear t1 

dt=0.01;		
% !计算时间（秒）
T=300;
NNT=T/dt;
t=0:dt:T;
P1=fangdaxishu*sin(THETA*t);
P2=fangdaxishu*sin(2*THETA*t)*50;
% degrees=200;
P_eachpoint=zeros(101,length(t));
P_eachpoint(50,:)=P1;
P_eachpoint(20,:)=P2;
P=P_eachpoint;


points=101;
PP=zeros(matrixsize,length(t));
for t1=1:matrixsize
        PP(t1,:)=Peq(t1,mode_vec,KMmapping,P_eachpoint,points,t);
end

u0=zeros(matrixsize,1);
udot0=zeros(matrixsize,1);
nbeta=0.25;
ngam=0.5;

[u1 udot1 u2dot1] = NewmarkInt(t,MM,CC,KK,PP,ngam,nbeta,u0,udot0);


pointnumber=30;%查看某个点的振动时程
phiResult=phiY(pointnumber,KMmapping,mode_vec,matrixsize);
Dis=zeros(1,length(t));

for t1=1:matrixsize
    Dis=Dis+phiResult(t1).*u1(t1,:);
end
figure()
plot(t,Dis,'b')
xlabel('Time [s]'); ylabel("displacement of point:"+num2str(pointnumber)); 
title("comparison of the displacement of midpoint calculated by matlab and ANSYS")
hold on 
plot(t(1:end/2),u(162,1:end/2),'r')




function result=Peq(Pmode,mode_vec,Mmapping,P_eachpoint,points,t)
    result=zeros(1,size(t,2));
    for t1=1:points
        if sum(and(Mmapping.Node==t1,Mmapping.DOF=='UY'))==0
            result=result+0*P_eachpoint(t1,:);
        else
            position_index=Mmapping.MatrixEqn(find(and(Mmapping.Node==t1,Mmapping.DOF=='UY')));
            result=result+mode_vec(position_index,Pmode)*P_eachpoint(t1,:);
        end
    end
end

function result=P_eq(mode,temp_vec,Matrix)
    vec=temp_vec(:,mode);
    result=vec'*Matrix*vec;
end

function result=phiY(node,Mmapping,mode_vec,nModes)
    position_index=Mmapping.MatrixEqn(find(and(Mmapping.Node==node,Mmapping.DOF=='UY')));
    result=mode_vec(position_index,1:nModes);
end
