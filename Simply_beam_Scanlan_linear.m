%% 导入数据
clc;clear all;close all
hb_to_mm ( 'KMatrix.matrix', 'K.txt' );
hb_to_mm ( 'MMatrix.matrix', 'M.txt' );

% scalevector=[1 10^-6];
scalevector=[1 10^-6];
scale2vector=[0.000 0.000];
Dis_result=[];
for k0=1:2
    scale=scalevector(k0);
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
% 还原对角线以上元素，使之为对称阵, ANSYS只给出下三角矩阵
K =diag(diag(Kmatrix)/2)+Kmatrix-diag(diag(Kmatrix));
K = K+K';
M =diag(diag(Mmatrix)/2)+Mmatrix-diag(diag(Mmatrix));
M = M+M';
% 特征值分析，即计算频率Freq和振型Phi，50代表求50阶，SM表示从较小的特征值开始求解
calmodes=50;
[eig_vec,eig_val]=eigs(K,M,calmodes,'SM');
[nfdof,nfdof]=size(eig_vec);
for j=1:nfdof
    mnorm=sqrt(eig_vec(:,j)'*M*eig_vec(:,j));
    eig_vec(:,j)=eig_vec(:,j)/mnorm;%振型质量归一化
end
[omeg,w_order] =sort(sqrt(diag(eig_val)));
mode_vec=eig_vec(:,w_order);
Freq = omeg/(2*pi);

save StiffMatrix.mat K M;
save mode_vec.mat mode_vec

%% Calculate the damping matrix of the beam
xi=0.003;
% If set alpha=0
alpha=0;
beta=2*xi./omeg;

C=zeros(calmodes,size(M,1),size(M,2));
for t1=1:calmodes
    C(t1,:,:)=alpha.*M+beta(t1).*K;
end

% beta=2*xi./omeg(1);
% beta=0;

% C=zeros(calmodes,size(M,1),size(M,2));
% for t1=1:calmodes
%     C(t1,:,:)=alpha.*M+beta.*K;
% end


%% TMD parameters
nTMD=2;
mTMD=[100 100]*scale;
% mTMD=[100 100]*10^-6;
cTMD=[2*mTMD(1)*omeg(1)*0.05 2*mTMD(2)*omeg(2)*0.05];
kTMD=[mTMD(1)*omeg(1)^2 mTMD(2)*omeg(1)^2];
nodeTMD=[50 20];   %Node number(location of the TMD)


%% Number of modes considered
nModes=50;
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
        C_temp=Cmatrix(t1,C,M);
        CC(t1,t1)=P_eq(t1,mode_vec,C_temp);
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

%% Scanlan linear model parameters
rho=1.225;
U=20.92*4.185*55/1000;
fn=4.185;
D=55/1000;
CL=0.613;
H1=-26.98;
omega=2*pi*fn;
Kw=omega*D/U;
C_aero=zeros(matrixsize,matrixsize);
for k1=1:nModes
    C_aero(k1,k1) = -0.5*rho*U*D*Kw*H1;
end
clear k1
CCa=CC+C_aero*scale2vector(k0);





%% Calculate the response
dt=0.01;
T=2500;
t=0:dt:T;
P0=0.5*rho*U^2*D*CL;
P0=P0*10000;
Pmode1=1;
points=101;
P_eachpoint1=P_mode_eachpoint(1,KMmapping,P0,points,omeg,t,mode_vec);
P_eachpoint=P_eachpoint1;

P=ones(matrixsize,length(t));

for t1=1:matrixsize
    if t1<=nModes
        P(t1,:)=Peq(t1,mode_vec,KMmapping,P_eachpoint,points,t);
    else
        P(t1,:)=zeros(1,size(t,2));
    end
end


nbeta=0.25;
ngam=0.5;
u0=zeros(matrixsize,1);
u0(1)=10000;
udot0=zeros(matrixsize,1);
[u udot u2dot] = NewmarkInt(t,MM,CC,KK,P,ngam,nbeta,u0,udot0);
[ua udota u2dota] = NewmarkInt(t,MM,CCa,KK,P,ngam,nbeta,u0,udot0);


% figure()
% grid on;
% plot(t,u(1,:),'b');
% xlabel('Time [s]'); ylabel('displacement of generalized coordinate'); 

% figure()
% plot(t,u(end-1,:),'b')
% xlabel('Time [s]'); ylabel('displacement of TMD'); 

% close all
% pointnumber=50;%查看某个点的振型
% phiResult=phiY(pointnumber,KMmapping,mode_vec,nModes);
% Dis=zeros(1,length(t));
% for t1=1:nModes
%     Dis=Dis+phiResult(t1).*u(t1,:);
% end
% 
% close all
pointnumber=50;%查看某个点的振型
phiResult=phiY(pointnumber,KMmapping,mode_vec,nModes);
Disa=zeros(1,length(t));
for t1=1:nModes
    Disa=Disa+phiResult(t1).*ua(t1,:);
end
% figure()
% plot(t,Dis,'b')
% xlabel('Time [s]'); ylabel("displacement of point:"+num2str(pointnumber)); 
% title("comparison of the displacement of midpoint calculated by matlab")
% hold on
% plot(t,Disa,'r')
% legend("without aerodynamic damping","with aerodynamic damping")

Dis_result=[Dis_result;Disa];
end
figure
plot(t,Dis_result(1,:),'b')
xlabel('Time [s]'); ylabel("displacement of point:"+num2str(pointnumber)); 
title("comparison of the displacement of midpoint calculated by matlab")
hold on
plot(t,Dis_result(2,:),'r')
% legend("aerodynamic damping with TMD","aerodynamic damping without TMD")
legend("no aerodynamic damping with TMD"," no aerodynamic damping without TMD")
% 
% [psd_avg, f, psd_plot]=fft_transfer(1/dt,Dis_result(1,:)');
% figure
% plot(f,psd_plot)

% figure
% plot(t,Dis_result(1,:),'b')
% xlabel('Time [s]'); ylabel("displacement of point:"+num2str(pointnumber)); 
% title("comparison of the displacement of midpoint calculated by matlab")
% hold on
% plot(t,Dis_result(2,:),'r')
% legend("TMD with aerodynamic damping","TMD without aerodynamic damping")
% 
% [psd_avg, f, psd_plot]=fft_transfer(1/dt,Dis_result(1,:)');
% figure
% plot(f,psd_plot)

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

function result=P_mode_eachpoint(Pmode,Mmapping,P0,points,omeg,t,mode_vec)
    for t1=1:points
        if sum(and(Mmapping.Node==t1,Mmapping.DOF=='UY'))==0
            result(t1,:)=0*P0*sin(omeg(Pmode)*t);
        else
            position_index=Mmapping.MatrixEqn(find(and(Mmapping.Node==t1,Mmapping.DOF=='UY')));
            result(t1,:)=mode_vec(position_index,Pmode)*P0*sin(omeg(Pmode)*t);
        end
    end
end



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

