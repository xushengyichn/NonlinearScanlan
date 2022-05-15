clc;clear all;close all
hb_to_mm ( 'KMatrix.matrix', 'K.txt' );
hb_to_mm ( 'MMatrix.matrix', 'M.txt' );
KMmapping = importmappingmatrix('KMatrix.mapping');

%%nodemumber:2001-2205
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
% 特征值分析，即计算频率Freq和振型Phi
[eig_vec,eig_val]=eigs(K,M,50,'SM');
[nfdof,nfdof]=size(eig_vec);
for j=1:nfdof
    mnorm=sqrt(eig_vec(:,j)'*M*eig_vec(:,j));
    eig_vec(:,j)=eig_vec(:,j)/mnorm;
end
[omeg,w_order] =sort(sqrt(diag(eig_val)));
mode_vec=eig_vec(:,w_order);
Freq = omeg/(2*pi);

save StiffMatrix.mat K M;