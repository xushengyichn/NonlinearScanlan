%% 导入数据
clc;clear all;close all
hb_to_mm ( 'KMatrix.matrix', 'K.txt' );
hb_to_mm ( 'MMatrix.matrix', 'M.txt' );

% scalevector=[1 10^-6];
scalevector=[1 10^-6];

Dis_result=[];
for k0=1:1
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
nTMD=0;
mTMD=[100]*scale;
% mTMD=[100 100]*10^-6;
cTMD=[2*mTMD(1)*omeg(1)*0.05];
kTMD=[mTMD(1)*omeg(1)^2];
nodeTMD=[50];   %Node number(location of the TMD)


%% Number of modes considered
nModes=1;
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

%% Scanlan nonlinear model parameters
rho=1.225;
U=20.92*4.185*55/1000;
fn=4.185;
D=55/1000;
omega=2*pi*fn;
Kw=omega*D/U;
Y1=0.1;
epsilon=30;
Y2=0;
gamma = 1/2; % Factor in the Newmark algorithm 
beta = 1/4; % Factor in the Newmark algorithm 





%% Calculate the response
h=0.01;
T=1000;
t=0:h:T;
P0=0.5*rho*U^2*D*0;
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
u0(1)=0.000001;
udot0=zeros(matrixsize,1);
pp=P;
gfun = @(u,udot) Scanlan_nonlinear(u,udot,MM,CC,KK,gamma,beta,h,rho,U,D,Y1,epsilon,Y2,matrixsize,nModes);
u = nonlinear_newmark_krenk(gfun,MM,pp,u0,udot0,gamma,beta,h);

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
Dis=zeros(1,length(t));
for t1=1:nModes
    Dis=Dis+phiResult(t1).*u(t1,:);
end
figure()
plot(t,Dis,'b')
xlabel('Time [s]'); ylabel("displacement of point:"+num2str(pointnumber)); 
title("comparison of the displacement of midpoint calculated by matlab")
% hold on
% plot(t,Disa,'r')
% legend("without aerodynamic damping","with aerodynamic damping")
close all
Dis_result=[Dis_result;Dis];
end
figure
plot(t,Dis_result(1,:),'b')
xlabel('Time [s]'); ylabel("displacement of point:"+num2str(pointnumber)); 
title("comparison of the displacement of midpoint calculated by matlab")
% hold on
% figure
% plot(t,u(end,:),'r')
% % legend("aerodynamic damping with TMD","aerodynamic damping without TMD")
% title("TMD")

% [psd_avg, f, psd_plot]=fft_transfer(1/h,Dis_result(1,:)');
% figure
% plot(f,psd_plot)
% xlim([0 0.3])




















%% below is the necessary functions
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

%% Nonlinear Newmark algorithm
function [u, udot, u2dot] = nonlinear_newmark_krenk(gfun,MM,pp,u0,udot0,gamma,beta,h)
    % Initialize variables
    u = zeros(size(MM,1),size(pp,2));
    udot = zeros(size(MM,1),size(pp,2));
    u2dot = zeros(size(MM,1),size(pp,2));
    % Initial conditions
    u(:,1) = u0; % Initial displacements
    udot(:,1) = udot0; % Initial velocities
    g = gfun(u0,udot0); % Evaluate the nonlinear function
    u2dot(:,1)=MM\(pp(:,1)-g); % Calculate the initial accelerations
    for ii=1:size(pp,2)-1
        % Prediction step
        u2dot(:,ii+1)=u2dot(:,ii); % Predicted accelerations
        udot(:,ii+1)=udot(:,ii)+h*u2dot(:,ii); % Predicted velocities 
        u(:,ii+1)=u(:,ii)+h*udot(:,ii)+1/2*h^2*u2dot(:,ii); % Predicted displacements
        Nit=0; % Number of iterations 
        konv=0; % Has the iterations converged 0=No, 1=Yes
        % Reusidal calculation, system matrices and increment correction
        while Nit<1000 && konv==0
            Nit=Nit+1;
            [g, Ks] = gfun(u(:,ii+1),udot(:,ii+1)); % Calculate function value and the tangent
            rr = pp(:,ii+1)-MM*u2dot(:,ii+1) - g; % Calculate residual
            du = Ks\rr; % Increment correction
            u(:,ii+1)=u(:,ii+1)+du; % Add incremental correction to the displacement
            udot(:,ii+1)=udot(:,ii+1)+gamma*h/(beta*h^2)*du; % Incremental correction for velocities 
            u2dot(:,ii+1)=u2dot(:,ii+1)+1/(beta*h^2)*du; % Incremental correction accelerations
            if sqrt(rr'*rr)/length(rr)<1.0e-8 % Convergence criteria
                konv=1; % konv = 1 if the convergence criteria is fulfilled.
            end
        end
    end
    end

    %% Function file for the Scanlan nonlinear model
function [g,ks]= Scanlan_nonlinear(u,udot,MM,CC,KK,gamma,beta,h,rho,U,D,Y1,epsilon,Y2,matrixsize,nModes)
    diagMatrix=zeros(matrixsize,matrixsize);
    for k1=1:nModes
        diagMatrix(k1,k1)=1;
    end

    g = (CC-diagMatrix*rho.*U.*D.*Y1.*(1-epsilon.*u.^2./D.^2))*udot+(KK-diagMatrix.*rho.*U.^2.*Y2./D)*u; % Function value
    kc= CC-diagMatrix*rho.*U.*D.*Y1.*(1-epsilon.*u.^2./D.^2);
    u_udot=zeros(matrixsize,matrixsize);
    for k1=1:size(u_udot,1)
        u_udot(k1,k1)=u(k1)*udot(k1);
    end
    kspring= KK-diagMatrix*rho.*U.^2.*Y2./D +2.*rho.*U.*D.*Y1.*epsilon*u_udot;
    ks = kspring + gamma.*h./(beta.*h.^2).*kc + 1./(beta.*h.^2).*MM; % Linearization
end