% clc
% clear
% close all
% 
% 
% % 调试参数
% % Debug parameter
% %% 将ANSYS中的稀疏矩阵处理为完全矩阵
% %% Handling sparse matrices in ANSYS as full matrices
% % 导入ANSYS MCK矩阵
% % Import MCK matrix from ANSYS
% hb_to_mm ('KMatrix.matrix', 'K.txt');
% hb_to_mm ('MMatrix.matrix', 'M.txt');
% hb_to_mm ('CMatrix.matrix', 'C.txt');
% 
% %map the node and matrix from the KMatrix.mapping and MMatrix.mapping
% Kdata = importdata('K.txt').data;
% Kmatrix = zeros(Kdata(1, 1), Kdata(1, 2));
% 
% for i = 2:size(Kdata, 1)
%     Kmatrix(Kdata(i, 1), Kdata(i, 2)) = Kdata(i, 3);
% end
% 
% Mdata = importdata('M.txt').data;
% Mmatrix = zeros(Mdata(1, 1), Mdata(1, 2));
% 
% for i = 2:size(Mdata, 1)
%     Mmatrix(Mdata(i, 1), Mdata(i, 2)) = Mdata(i, 3);
% end
% 
% Cdata = importdata('C.txt').data;
% Cmatrix_DP = zeros(Cdata(1, 1), Cdata(1, 2));
% 
% for i = 2:size(Cdata, 1)
%     Cmatrix_DP(Cdata(i, 1), Cdata(i, 2)) = Cdata(i, 3);
% end
% 
% % 还原对角线以上元素，使之为对称阵, ANSYS只给出下三角矩阵
% % Restore the elements above the diagonal to make it a symmetric matrix, ANSYS only gives the lower triangular matrix
% K = diag(diag(Kmatrix) / 2) + Kmatrix - diag(diag(Kmatrix));
% K = K + K';
% M = diag(diag(Mmatrix) / 2) + Mmatrix - diag(diag(Mmatrix));
% M = M + M';
% C_exp = diag(diag(Cmatrix_DP) / 2) + Cmatrix_DP - diag(diag(Cmatrix_DP));
% C_exp = C_exp + C_exp';
% C =C_exp;
% 
% 
% 
% % 设置模态叠加法生成阻尼矩阵的方式
% % 1 表示采用瑞利阻尼矩阵生成对应的模态叠加法阻尼矩阵
% % 2 表示采用模态阻尼比直接生成的模态叠加法阻尼矩阵，需要指定模态阻尼比zeta
% % Set the method for generating the damping matrix by the modal superposition method
% % 1 means use the Rayleigh damping matrix to generate the corresponding modal superposition method damping matrix
% % 2 means use the modal damping ratio to directly generate the modal superposition method damping matrix, you need to specify the modal Damping ratio zeta
% DampingMatrixParameter = 2;
% 
% % 设置桥面坐标导入模式
% % 1 采用程序自动写入一个含有桥面节点的数组，从 1 到 n
% % 2 读入文件
% % Set the bridge deck coordinate import mode
% % 1 Use the program to automatically write an array containing the bridge deck nodes, from 1 to n
% % 2 read in file
% nodeondeckimport = 1;
% % 设置外荷载施加方式
% % 1 施加模态力
% % 2 施加节点力（需要调整对应施加荷载部分的代码，未完成）
% % Set the external load application method
% % 1 Apply modal force
% % 2 Apply nodal forces (need to adjust the code corresponding to the applied load section)
% externalforcemethod = 1;
% % 是否输出振动视频
% % Whether to output vibration video
% % 1 export vides
% % 2 do not export vides
% exportvideo=2;
% 
% % Number of modes considered
% nModes = 10; %考虑模态
% 
% dt = 0.1;
% % 计算时间（秒）Computation time (seconds)
% T_default = 2000;
% 
% m_range=linspace(100,100,1);
% zeta_range=linspace(0.01,0.1,10);
% f_range=linspace(0.05,0.6500,100);
% position_range=linspace(2,100,50);
% p_range=linspace(1,3,3);
% 
% Cases = combvec(m_range,zeta_range,f_range,position_range,p_range)';
% 
% omega_range=2*pi*Cases(:,3);
% k_range=Cases(:,1).*omega_range.^2;
% c_range=2.*Cases(:,1).*omega_range.*Cases(:,2);
% 
% TestingCases=[Cases(:,1) c_range k_range Cases(:,4) Cases(:,5)];
% 
% 
% ReferenceData = [];
% 
% for k1=1:size(p_range,2)
%     nTMD = 0;
%     mTMD = 0;
%     cTMD = 0;
%     kTMD = 0;
%     nodeTMD = 0; %Node number(location of the TMD)
%     MFC=[zeros(k1-1,1);1];
%     T =T_default;
%     [Dis,flag]=NonlinearScanlan_Beam_function(K,C,M,nTMD, mTMD, cTMD, ...
%     kTMD, nodeTMD, nModes,DampingMatrixParameter,nodeondeckimport, ...
%     externalforcemethod,exportvideo,dt,T,"zeta",0.3/100,"MFC",MFC,"points",101);
% 
%     while flag ==1
%         T=T+2000;
%         disp("计算未收敛，调整计算时间长度为"+num2str(T)+"秒")
%         [Dis,flag]=NonlinearScanlan_Beam_function(K,C,M,nTMD, mTMD, cTMD, ...
%     kTMD, nodeTMD, nModes,DampingMatrixParameter,nodeondeckimport, ...
%     externalforcemethod,exportvideo,dt,T,"zeta",0.3/100,"MFC",MFC,"points",101);
%     end
%     
% 
%     Dis_sel=Dis(:,end-round(end/10,0):end);
%     Max_Value=max(max(abs(Dis_sel)));
%     [Max_Position,Max_time]=find(abs(Dis_sel)==Max_Value);
% 
%     ReferenceData(k1,:)=[k1 Max_Position Max_Value];
%     clear Max_Position Max_time Max_Value
% end
% 
% 
% tic
% 
% collectdata1=zeros(size(TestingCases,1),1);
% collectdata2=zeros(size(TestingCases,1),1);
% collectdata3=zeros(size(TestingCases,1),1);
% collectdata4=zeros(size(TestingCases,1),1);
% collectdata5=zeros(size(TestingCases,1),1);
% 
% numIterations=100;    %total number of parfor iterations
% ppm = ParforProgressbar(numIterations);
% 
% parfor k1=1:100
% % parfor k1=1:size(TestingCases,1)
%     nTMD = 1;
%     mTMD = TestingCases(k1,1);
%     cTMD = TestingCases(k1,2);
%     kTMD = TestingCases(k1,3);
%     nodeTMD = TestingCases(k1,4); %Node number(location of the TMD)
%     MFC=[zeros(TestingCases(k1,5)-1,1);1];
%     T =T_default;
%     [Dis,flag]=NonlinearScanlan_Beam_function(K,C,M,nTMD, mTMD, cTMD, ...
%     kTMD, nodeTMD, nModes,DampingMatrixParameter,nodeondeckimport, ...
%     externalforcemethod,exportvideo,dt,T,"zeta",0.3/100,"MFC",MFC,"points",101);
% 
%     while flag ==1
%         T=T+2000;
%         disp("计算未收敛，调整计算时间长度为"+num2str(T)+"秒")
%         [Dis,flag]=NonlinearScanlan_Beam_function(K,C,M,nTMD, mTMD, cTMD, ...
%     kTMD, nodeTMD, nModes,DampingMatrixParameter,nodeondeckimport, ...
%     externalforcemethod,exportvideo,dt,T,"zeta",0.3/100,"MFC",MFC,"points",101);
%     end
%     
%     maxdis=max(max(Dis));
%     [Max_Position,Max_time]=find(Dis==maxdis);
%     maxref=max(Dis(ReferenceData(TestingCases(k1,5),2),:));
%     reference=TestingCases(k1,5);
% %     collectdata(k1,1:5)=[ReferenceData(reference,2) ReferenceData(reference,3) maxref Max_Position maxdis];
%     collectdata1(k1,1)=[ReferenceData(reference,2)];
%     collectdata2(k1,1)=[ReferenceData(reference,3)];
%     collectdata3(k1,1)=[maxref];
%     collectdata4(k1,1)=[Max_Position];
%     collectdata5(k1,1)=[maxdis];
% 
%        ppm.increment();
% 
% end
% toc
% delete(ppm);
% TestingCases = [TestingCases collectdata1 collectdata2 collectdata3 collectdata4  collectdata5];
% Title = ["m","c","k","position","pmode","Original_max_point",...
%     "Corresponding_Dis","Dis_with_TMD","New_max_point","New_Dis"];
% Results=array2table(TestingCases,"VariableNames",Title) ;
% Results.Improvement=(Results.Corresponding_Dis-Results.Dis_with_TMD)./Results.Corresponding_Dis*100;
% 
% % save result Results

clc
clear
close all

load result
Results.Diff=(Results.New_Dis-Results.Dis_with_TMD)./Results.Dis_with_TMD*100;
m=100;
data=[Results.k,Results.c,Results.position,Results.pmode,Results.Improvement];
f=sqrt(data(:,1)./m)/2/pi;
omega= 2*pi*f;
zeta= data(:,2)./(2*m.*omega);
data2=[f,zeta, Results.position,Results.pmode,Results.Improvement];



% 点的大小是CPEE*100，对CPEE着色，'.'表示点的形状
mode=3;
data3=data2(data2(:,4)==mode,:);

scatter3(data3(:,1),data3(:,2),data3(:,3),data3(:,4)*10,data3(:,5),'.')
xlabel('Frequency')
ylabel('Zeta')
zlabel('Position')
grid on
h = colorbar;% 右侧颜色栏
set(get(h,'label'),'string','Improvement');% 给右侧颜色栏命名
title("Control effects with third-order modal aerodynamic force")
% xlim([0 20]) % X,Y,Z轴取值范围
% ylim([0 0.011])
% zlim([0 4])


% plot(Results.Diff)
% title("Max accleration vs Acceleration at the original point with maximum acceleration")
% xlabel("Testing cases")
% ylabel("Increase factor")

