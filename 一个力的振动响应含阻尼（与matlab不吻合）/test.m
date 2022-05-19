clc
clear
close all

load vectest
KMmapping = importmappingmatrix('KMatrix.mapping');
KMmapping3 = importmappingmatrix('KMatrixTMD.mapping');

TMDnum=[10001 10002];
V4=[];

for k1=1:size(KMmapping3,1)
    if k1<=size(KMmapping,1)
        Nodetemp=KMmapping.Node(k1);
        DOFtemp=KMmapping.DOF(k1);
        position_index=KMmapping3.MatrixEqn(find(and(KMmapping3.Node== Nodetemp,KMmapping3.DOF==DOFtemp)));
        V4(k1,:)=V3(position_index,:);
    else
        Nodetemp=TMDnum(KMmapping.Node(k1-size(KMmapping,1)));
        DOFtemp='UY';
        position_index=KMmapping3.MatrixEqn(find(and(KMmapping3.Node== Nodetemp,KMmapping3.DOF==DOFtemp)));
        V4(k1,:)=V3(position_index,:);
    end
end
clear k1

for k1=1:size(V2,2)
    if V2(end,k1)/V4(end,k1)<0
        V2(:,k1)=V2(:,k1)*-1;
    end
end

V1new=V1(1:end,1:(size(V2,2)-size(TMDnum,2)));
V2new=V2(1:end-2,:);
V3new=V4(1:end-2,1:(size(V2,2)-size(TMDnum,2)));

clearvars -except V1new V2new V3new
V1_line=90;
matrix=zeros(size(V1new,2));
list=zeros(size(V1new,2),1);
for k1 = 1:size(matrix,1)
    for k2 = 1:size(matrix,2)
        matrix(k1,k2)=V1new(V1_line,k2)*V2new(k2,k1);
    end
end
clear k1 k2
for k1 = 1:size(matrix,1)
    list(k1)=V3new(V1_line,k1);
end

x=matrix^(-1)*list

