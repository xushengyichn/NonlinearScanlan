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





clearvars -except K M C


hb_to_mm ( 'KMatrixTMD1.matrix', 'KTMD1.txt' );
hb_to_mm ( 'MMatrixTMD1.matrix', 'MTMD1.txt' );
hb_to_mm ( 'CMatrixTMD1.matrix', 'CTMD1.txt' );



Kdata = importdata('KTMD1.txt').data;
Kmatrix = zeros(Kdata(1,1),Kdata(1,2));
for i = 2:size(Kdata,1)
    Kmatrix(Kdata(i,1),Kdata(i,2)) = Kdata(i,3);
end

Mdata = importdata('MTMD1.txt').data;
Mmatrix = zeros(Mdata(1,1),Mdata(1,2));

for i = 2:size(Mdata,1)
    Mmatrix(Mdata(i,1),Mdata(i,2)) = Mdata(i,3);
end

Cdata = importdata('CTMD1.txt').data;
Cmatrix = zeros(Cdata(1,1),Cdata(1,2));

for i = 2:size(Cdata,1)
    Cmatrix(Cdata(i,1),Cdata(i,2)) = Cdata(i,3);
end

% 还原对角线以上元素，使之为对称阵, ANSYS只给出下三角矩阵
K1 =diag(diag(Kmatrix)/2)+Kmatrix-diag(diag(Kmatrix));
K1 = K1+K1';
M1 =diag(diag(Mmatrix)/2)+Mmatrix-diag(diag(Mmatrix));
M1 = M1+M1';
C1 =diag(diag(Cmatrix)/2)+Cmatrix-diag(diag(Cmatrix));
C1 = C1+C1';


KMmapping = importmappingmatrix('KMatrixTMD.mapping');
KMmapping1 = importmappingmatrix('KMatrixTMD1.mapping');























% 
% 
% 
% fangdaxishu=10;
% PI=3.14159265359;
% FRE=0.069704171453635;	
% 
% THETA=2*PI*FRE;
% 
% dt=0.01;		
% % !计算时间（秒）
% T=300;
% NNT=T/dt;
% t=0:dt:T;
% P1=fangdaxishu*sin(THETA*t);
% degrees=size(M,1);
% P_eachpoint=zeros(degrees,length(t));
% P_eachpoint(1,:)=P1;
% P=P_eachpoint;
% % P=zeros(50,length(t));
% % for t1=1:matrixsize
% %     if t1<=nModes
% %         P(t1,:)=Peq(t1,mode_vec3,KMmapping,P_eachpoint,points,t);
% %     else
% %         P(t1,:)=zeros(1,size(t,2));
% %     end
% % end
% KMmapping = importmappingmatrix('KMatrixTMD.mapping');
% u0=zeros(degrees,1);
% udot0=zeros(degrees,1);
% nbeta=0.25;
% ngam=0.5;
% 
% [u udot u2dot] = NewmarkInt(t,M,C,K,P,ngam,nbeta,u0,udot0);
% clear u0 udot0 P_eachpoint P1
% 
% figure()
% plot(t,u(1,:),'r')
% xlabel('Time [s]'); ylabel("displacement of point"); 
% title("comparison of the displacement of midpoint calculated by matlab and ANSYS")
% hold on 
% 
% 
% dataANSYS=readmatrix("T_DIS50.txt");
% dataANSYS=dataANSYS(2:end,:);
% plot(dataANSYS(:,1),dataANSYS(:,2),'b')
% 
% 
% figure()
% plot(t,u(141,:),'r')
% xlabel('Time [s]'); ylabel("displacement of point"); 
% title("comparison of the displacement of midpoint calculated by matlab and ANSYS")
% hold on 
% 
% 
% dataANSYS=readmatrix("T_DIS20.txt");
% dataANSYS=dataANSYS(2:end,:);
% plot(dataANSYS(:,1),dataANSYS(:,2),'b')
% 
% figure()
% plot(t,u(143,:),'r')
% xlabel('Time [s]'); ylabel("displacement of point"); 
% title("comparison of the displacement of midpoint calculated by matlab and ANSYS")
% hold on 
% 
% 
% dataANSYS=readmatrix("TMD2_DIS.txt");
% dataANSYS=dataANSYS(2:end,:);
% plot(dataANSYS(:,1),dataANSYS(:,2),'b')
% 
% 
% figure()
% plot(t,u(202,:),'r')
% xlabel('Time [s]'); ylabel("displacement of point"); 
% title("comparison of the displacement of midpoint calculated by matlab and ANSYS")
% hold on 
% 
% 
% dataANSYS=readmatrix("TMD1_DIS.txt");
% dataANSYS=dataANSYS(2:end,:);
% plot(dataANSYS(:,1),dataANSYS(:,2),'b')



