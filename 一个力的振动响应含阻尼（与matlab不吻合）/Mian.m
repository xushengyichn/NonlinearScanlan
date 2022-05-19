%% 导入数据
clc;clear;close all
hb_to_mm ( 'KMatrix.matrix', 'K.txt' );
hb_to_mm ( 'MMatrix.matrix', 'M.txt' );
hb_to_mm ( 'CMatrix.matrix', 'C.txt' );


%%map the node and matrix from the KMatrix.mapping and MMatrix.mapping 
Kdata = importdata('K.txt').data;
Kmatrix = zeros(Kdata(1,1),Kdata(1,2));
for i = 2:size(Kdata,1)
    Kmatrix(Kdata(i,1),Kdata(i,2)) = Kdata(i,3);
end

Mdata = importdata('M.txt').data;
Mmatrix = zeros(Mdata(1,1),Mdata(1,2));

for i = 2:size(Mdata,1)
    Mmatrix(Mdata(i,1),Mdata(i,2)) = Mdata(i,3);
end

Cdata = importdata('C.txt').data;
Cmatrix_DP = zeros(Cdata(1,1),Cdata(1,2));

for i = 2:size(Cdata,1)
    Cmatrix_DP(Cdata(i,1),Cdata(i,2)) = Cdata(i,3);
end
% 还原对角线以上元素，使之为对称阵, ANSYS只给出下三角矩阵
K =diag(diag(Kmatrix)/2)+Kmatrix-diag(diag(Kmatrix));
K = K+K';
M =diag(diag(Mmatrix)/2)+Mmatrix-diag(diag(Mmatrix));
M = M+M';
C_exp =diag(diag(Cmatrix_DP)/2)+Cmatrix_DP-diag(diag(Cmatrix_DP));
C_exp = C_exp+C_exp';
C_exp = 0.01*K;
% 特征值分析，即计算频率Freq和振型Phi，50代表求50阶，SM表示从较小的特征值开始求解
calmodes=200;%考虑模态
[eig_vec,eig_val]=eigs(K,M,calmodes,'SM');
[nfdof,nfdof]=size(eig_vec);
for j=1:nfdof
    mnorm=sqrt(eig_vec(:,j)'*M*eig_vec(:,j));
    eig_vec(:,j)=eig_vec(:,j)/mnorm;%振型质量归一化
end
[omeg,w_order] =sort(sqrt(diag(eig_val)));
mode_vec=eig_vec(:,w_order);
Freq = omeg/(2*pi);

% % save StiffMatrix.mat K M;
% % save mode_vec.mat mode_vec

KMmapping = importmappingmatrix('KMatrix.mapping');



clearvars -except M C K Freq omeg mode_vec C_exp KMmapping

%% TMD parameters
nTMD=2;
mTMD=[100 100];
% mTMD=[100 100]*10^-6;
cTMD=[2*mTMD(1)*omeg(1)*0.05 2*mTMD(2)*omeg(2)*0.05];
kTMD=[mTMD(1)*omeg(1)^2 mTMD(2)*omeg(2)^2];
nodeTMD=[50 20];   %Node number(location of the TMD)


%% Number of modes considered
nModes=10;%考虑模态
matrixsize=nTMD+nModes;
KMmapping = importmappingmatrix('KMatrix.mapping');

phiTMD=zeros(nTMD,nModes);
% phiTMD row:TMD for each loaction column:the mode shape at the each
% location of tmd
for t1=1:nTMD
    for t2=1:nModes
        position_index=KMmapping.MatrixEqn(find(and(KMmapping.Node==nodeTMD(t1),KMmapping.DOF=='UY')));
        phiTMD(t1,t2)=mode_vec(position_index,t2);
    end
end
clear t1 t2
%% Create the Mass matrix
MM=zeros(matrixsize,matrixsize);
for t1=1:matrixsize
    if t1<=nModes
        MM(t1,t1)=P_eq(t1,mode_vec,M);
    end
    if and(t1>nModes,t1<=matrixsize)
        MM(t1,t1)=mTMD(t1-nModes);
    end
end
clear t1 

%% Create the Damping matrix

CC=zeros(matrixsize,matrixsize);
for t1=1:matrixsize
    if t1<=nModes
        CC(t1,t1)=P_eq(t1,mode_vec,C_exp);%P=parameters
    end
    if and(t1>nModes,t1<=matrixsize)
        CC(t1,t1)=cTMD(t1-nModes);
    end
end
clear t1


for t1=1:nModes
    for t2=1:nTMD
        for t3=1:nModes
            temp_data=cTMD(t2)*phiTMD(t2,t3)*phiTMD(t2,t1);
            CC(t1,t3)=CC(t1,t3)+temp_data;
            clear temp_data
        end
        indextmdt2=nModes+t2;
        CC(t1,indextmdt2)=CC(t1,indextmdt2)-cTMD(t2)*phiTMD(t2,t1);
    end
end
clear t1 t2 t3

for t1=1:nTMD
    for t2=1:nModes
        temp_data=-cTMD(t1)*phiTMD(t1,t2);
        CC(nModes+t1,t2)=CC(nModes+t1,t2)+temp_data;
        clear temp_data
    end
end
clear t1


% 
% CC=zeros(matrixsize,matrixsize);
% for t1=1:matrixsize
%     if t1<=nModes
%         C_temp=Cmatrix(t1,C,M);
%         CC(t1,t1)=P_eq(t1,mode_vec,C_temp);
%     end
%     if and(t1>nModes,t1<=matrixsize)
%         CC(t1,t1)=cTMD(t1-nModes);
%     end
% end
% clear t1 
% 
% 
% for t1=1:nModes
%     for t2=1:nTMD
%         for t3=1:nModes
%             temp_data=cTMD(t2)*phiTMD(t2,t3)*phiTMD(t2,t1);
%             CC(t1,t3)=CC(t1,t3)+temp_data;
%             clear temp_data
%         end
%         indextmdt2=nModes+t2;
%         CC(t1,indextmdt2)=CC(t1,indextmdt2)-cTMD(t2)*phiTMD(t2,t1);
%     end
% end
% clear t1 t2 t3
% 
% for t1=1:nTMD
%     for t2=1:nModes
%         temp_data=-cTMD(t1)*phiTMD(t1,t2);
%         CC(nModes+t1,t2)=CC(nModes+t1,t2)+temp_data;
%         clear temp_data
%     end
% end
% clear t1



%% Create the Stiffness Matrix

KK=zeros(matrixsize,matrixsize);
for t1=1:matrixsize
    if t1<=nModes
        KK(t1,t1)=P_eq(t1,mode_vec,K);%P=parameters
    end
    if and(t1>nModes,t1<=matrixsize)
        KK(t1,t1)=kTMD(t1-nModes);
    end
end
clear t1


for t1=1:nModes
    for t2=1:nTMD
        for t3=1:nModes
            temp_data=kTMD(t2)*phiTMD(t2,t3)*phiTMD(t2,t1);
            KK(t1,t3)=KK(t1,t3)+temp_data;
            clear temp_data
        end
        indextmdt2=nModes+t2;
        KK(t1,indextmdt2)=KK(t1,indextmdt2)-kTMD(t2)*phiTMD(t2,t1);
    end
end
clear t1 t2 t3

for t1=1:nTMD
    for t2=1:nModes
        temp_data=-kTMD(t1)*phiTMD(t1,t2);
        KK(nModes+t1,t2)=KK(nModes+t1,t2)+temp_data;
        clear temp_data
    end
end
clear t1


%% 特征值分析
clearvars -except KK MM CC matrixsize nModes nTMD mode_vec KMmapping
calmodes=matrixsize;%考虑模态
[eig_vec,eig_val]=eigs(KK,MM,calmodes,'SM');
[nfdof,nfdof]=size(eig_vec);
for k1=1:nfdof
    mnorm=sqrt(eig_vec(:,k1)'*MM*eig_vec(:,k1));
%     disp(mnorm)
    eig_vec(:,k1)=eig_vec(:,k1)/mnorm;%振型质量归一化
end
clear k1 ndof
[omeg,w_order] =sort(sqrt(diag(eig_val)));
mode_vec2=eig_vec(:,w_order);
Freq2 = omeg/(2*pi);

%% 导入ansys增加TMD后的矩阵

hb_to_mm ( 'KMatrixTMD.matrix', 'KTMD.txt' );
hb_to_mm ( 'KMatrixTMDre.matrix', 'KTMD2.txt' );
hb_to_mm ( 'MMatrixTMD.matrix', 'MTMD.txt' );
hb_to_mm ( 'CMatrixTMD.matrix', 'CTMD.txt' );

Kdata2 = importdata('KTMD2.txt').data;
Kmatrix2 = zeros(Kdata2(1,1),Kdata2(1,2));
for i = 2:size(Kdata2,1)
    Kmatrix2(Kdata2(i,1),Kdata2(i,2)) = Kdata2(i,3);
end
K2 =diag(diag(Kmatrix2)/2)+Kmatrix2-diag(diag(Kmatrix2));
K2 = K2+K2';

%%map the node and matrix from the KMatrix.mapping and MMatrix.mapping 
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
Cmatrix_DP = zeros(Cdata(1,1),Cdata(1,2));

for i = 2:size(Cdata,1)
    Cmatrix_DP(Cdata(i,1),Cdata(i,2)) = Cdata(i,3);
end
% 还原对角线以上元素，使之为对称阵, ANSYS只给出下三角矩阵
K =diag(diag(Kmatrix)/2)+Kmatrix-diag(diag(Kmatrix));
K = K+K';
M =diag(diag(Mmatrix)/2)+Mmatrix-diag(diag(Mmatrix));
M = M+M';
C =diag(diag(Cmatrix_DP)/2)+Cmatrix_DP-diag(diag(Cmatrix_DP));
C = C+C'+0.01*K2;

% 特征值分析，即计算频率Freq和振型Phi，50代表求50阶，SM表示从较小的特征值开始求解
calmodes=matrixsize;%考虑模态
[eig_vec,eig_val]=eigs(K,M,calmodes,'SM');
[nfdof,nfdof]=size(eig_vec);
for j=1:nfdof
    mnorm=sqrt(eig_vec(:,j)'*M*eig_vec(:,j));
    eig_vec(:,j)=eig_vec(:,j)/mnorm;%振型质量归一化
end
[omeg,w_order] =sort(sqrt(diag(eig_val)));
mode_vec3=eig_vec(:,w_order);
Freq3 = omeg/(2*pi);



matrixsize=nTMD+nModes;
KMmapping3 = importmappingmatrix('KMatrixTMD.mapping');
phiTMD3=zeros(nTMD,matrixsize);
nodeTMD3=[10001 10002];
% phiTMD row:TMD for each loaction column:the mode shape at the each
% location of tmd
for t1=1:nTMD
    for t2=1:matrixsize
        position_index=KMmapping3.MatrixEqn(find(and(KMmapping3.Node==nodeTMD3(t1),KMmapping3.DOF=='UY')));
        phiTMD3(t1,t2)=mode_vec3(position_index,t2);
    end
end
clear t1 t2

V1=mode_vec;
V2=mode_vec2;
V3=mode_vec3;
save vectest V1 V2 V3

%调整ANSYS振型使之与matlab对应
decidevalue2=zeros(size(phiTMD3,2),1);
for t1=1:size(phiTMD3,2)
    decidevalue=phiTMD3(1,t1)/mode_vec2(end-1,t1);
    temp=abs(mode_vec2(end-1,t1))-abs(phiTMD3(1,t1));
    if temp>10e-6
        disp("模态"+num2str(t1)+"的完全矩阵和缩减矩阵tmd振型差异过大"+num2str(temp))
    end
    if decidevalue<0
%         phiTMD3(:,t1)=phiTMD3(:,t1)*-1;
%         mode_vec3(:,t1)=mode_vec3(:,t1)*-1;
        mode_vec2(:,t1)=mode_vec2(:,t1)*-1;
        decidevalue2(t1)=1;
    end
end



%% Calculate the response
points=101;
fangdaxishu=10;
fangdaxishu2=100;
zhenxing1=0.99950656036586;
zhenxing2=-0.92977648588903;
PI=3.14159265359;
FRE=0.069704171453635;	
% FRE=Freq2(1);	
FRE2=0.278900577315756;
THETA=2*PI*FRE;
THETA2=2*PI*FRE2;
dt=0.1;		
% !计算时间（秒）
T=500;
NNT=T/dt;
t=0:dt:T;
P1=fangdaxishu*zhenxing1*sin(THETA*t);

P2=fangdaxishu2*zhenxing2*sin(THETA2*t)*10;
P_eachpoint=zeros(points,length(t));
P_eachpoint(50,:)=P1;
P_eachpoint(20,:)=P2;
PP=zeros(matrixsize,length(t));

for t1=1:matrixsize
    if t1<=nModes
        PP(t1,:)=Peq(t1,mode_vec,KMmapping,P_eachpoint,points,t);
    else
        PP(t1,:)=zeros(1,size(t,2));
    end
end


nbeta=0.25;
ngam=0.5;
u0=zeros(matrixsize,1);
udot0=zeros(matrixsize,1);
tic
[u1 udot u2dot] = NewmarkInt(t,MM,CC,KK,PP,ngam,nbeta,u0,udot0);
toc
% u1(decidevalue2==1,:)=u1(decidevalue2==1,:)*-1;%使位移与振型相对应
% pointnumber=50;%查看某个点的振动时程
% phiResult=phiY(pointnumber,KMmapping3,mode_vec3,nModes);
% for k1 = 1:length(decidevalue2)
%     if decidevalue2(k1)==1
%         phiResult(k1)=phiResult(k1)*-1;
%     end
% end
% Dis=zeros(1,length(t));
% for t1=1:nModes
%     Dis=Dis+phiResult(t1).*u1(t1,:);
% end

u0=zeros(202,1);
udot0=zeros(202,1);
P = zeros(202,length(t));
P(38,:) = P2;
P(99,:) = P1;
tic
[u udot u2dot] = NewmarkInt(t,M,C,K,P,ngam,nbeta,u0,udot0);
toc

% % figure()
% plot(t,Dis,'b')
% xlabel('Time [s]'); ylabel("displacement of point:"+num2str(pointnumber)); 
% title("comparison of the displacement of midpoint calculated by matlab and ANSYS")
% hold on
% plot(t,u(1,:),'r')
% legend("缩减","完整")


close all
% % 若只在tmd上作用荷载
% P_eachpoint=zeros(points,length(t));
% 
% PP=zeros(matrixsize,length(t));
% PP(end-1,:)=P1+P2;
% u0(end-1)=0.1;
% tic
% [u1 udot u2dot] = NewmarkInt(t,MM,CC,KK,PP,ngam,nbeta,u0,udot0);
% toc
% u1(decidevalue2==1,:)=u1(decidevalue2==1,:)*-1;%使位移与振型相对应
% 
% 
% u0=zeros(202,1);
% u0(101)=0.1;
% udot0=zeros(202,1);
% P = zeros(202,length(t));
% P(101,:) = P1+P2;
% tic
% [u udot u2dot] = NewmarkInt(t,M,C,K,P,ngam,nbeta,u0,udot0);
% toc
figure
plot(t,u1(end-1,:),'r')
hold on
plot(t(1:end*0.9),u(101,1:end*0.9),'b')
xlabel('Time [s]'); ylabel("displacement"); 
title("comparison of the displacement of TMD1 calculated by FUll matrix or superposition method")

figure
plot(t,u1(end,:),'r')
hold on
plot(t(1:end*0.9),u(40,1:end*0.9),'b')
xlabel('Time [s]'); ylabel("displacement"); 
title("comparison of the displacement of TMD2 calculated by FUll matrix or superposition method")

% figure()
% grid on;
% plot(t,u(1,:),'b');
% xlabel('Time [s]'); ylabel('displacement of generalized coordinate'); 
% 
% figure()
% plot(t,u(end-1,:),'b')
% xlabel('Time [s]'); ylabel('displacement of TMD'); 

% close all


% 
% pointnumber=50;%查看某个点的振动时程
% phiResult=phiY(pointnumber,KMmapping,mode_vec3,nModes);
% Dis=zeros(1,length(t));
% for t1=1:nModes
%     Dis=Dis+phiResult(t1).*u(t1,:);
% end
% % figure()
% plot(t,Dis,'b')
% xlabel('Time [s]'); ylabel("displacement of point:"+num2str(pointnumber)); 
% title("comparison of the displacement of midpoint calculated by matlab and ANSYS")
% hold on
% dataANSYS=readmatrix("T_DIS50.txt");
% dataANSYS=dataANSYS(2:end,:);
% plot(dataANSYS(:,1),dataANSYS(:,2),'r')
% 
% close all
pointnumber=50;%查看某个点的振动时程
phiResult=phiY(pointnumber,KMmapping,mode_vec,nModes);
Dis=zeros(1,length(t));
for t1=1:nModes
    Dis=Dis+phiResult(t1).*u1(t1,:);
end
figure()
plot(t,Dis,'b')
xlabel('Time [s]'); ylabel("displacement of point:"+num2str(pointnumber)); 
title("comparison of the displacement calculated by FUll matrix or superposition method")
hold on
plot(t(1:end*0.9),u(99,1:end*0.9),'r')


pointnumber=30;%查看某个点的振动时程
phiResult=phiY(pointnumber,KMmapping,mode_vec,nModes);
Dis=zeros(1,length(t));
for t1=1:nModes
    Dis=Dis+phiResult(t1).*u1(t1,:);
end
figure()
plot(t,Dis,'b')
xlabel('Time [s]'); ylabel("displacement of point:"+num2str(pointnumber)); 
title("comparison of the displacement calculated by FUll matrix or superposition method")
hold on
plot(t(1:end*0.9),u(59,1:end*0.9),'r')

% 
% close all
% 
% pointnumber=10001;%查看某个点的振动时程
% 
% figure()
% plot(t,u(end-1,:),'b')
% xlabel('Time [s]'); ylabel("displacement of point:"+num2str(pointnumber)); 
% title("comparison of the displacement of TMD1")
% hold on
% dataANSYS=readmatrix("TMD1_DIS.txt");
% dataANSYS=dataANSYS(2:end,:);
% plot(dataANSYS(:,1),dataANSYS(:,2),'r')
% close all
% 
% pointnumber=10002;%查看某个点的振动时程
% figure()
% plot(t,u(end,:),'b')
% xlabel('Time [s]'); ylabel("displacement of point:"+num2str(pointnumber)); 
% title("comparison of the displacement of TMD2")
% hold on
% dataANSYS=readmatrix("TMD2_DIS.txt");
% dataANSYS=dataANSYS(2:end,:);
% plot(dataANSYS(:,1),dataANSYS(:,2),'r')
% 
% close all










function result=P_eq(mode,temp_vec,Matrix)
    vec=temp_vec(:,mode);
    result=vec'*Matrix*vec;
end

function result=Cmatrix(mode,C,referencematrix)
    size_temp=size(referencematrix,1);
    C_temp=C(mode,:,:);
    result=reshape(C_temp,[size_temp,size_temp]);
end

function result=phiY(node,Mmapping,mode_vec,nModes)
    position_index=Mmapping.MatrixEqn(find(and(Mmapping.Node==node,Mmapping.DOF=='UY')));
    result=mode_vec(position_index,1:nModes);
end

% function result=P_mode_eachpoint(Pmode,Mmapping,P0,points,omeg,t,mode_vec)
%     for t1=1:points
%         if sum(and(Mmapping.Node==t1,Mmapping.DOF=='UY'))==0
%             result(t1,:)=0*P0*sin(omeg(Pmode)*t);
%         else
%             position_index=Mmapping.MatrixEqn(find(and(Mmapping.Node==t1,Mmapping.DOF=='UY')));
%             result(t1,:)=mode_vec(position_index,Pmode)*P0*sin(omeg(Pmode)*t);
%         end
%     end
% end

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
