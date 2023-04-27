%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%Date: 2022-09-12 10:58:11
%LastEditors: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%LastEditTime: 2022-09-12 11:20:42
%FilePath: \NonlinearScanlan\CreateMatrixwithTMD.m
%Description: 生成节段模型，安装TMD后的质量、阻尼和刚度矩阵(不考虑气动阻尼和气动刚度)
%
%Copyright (c) 2022 by xushengyichn 54436848+xushengyichn@users.noreply.github.com, All Rights Reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [MM, CC, KK] = CreateMatrixwithTMD(nModes,Mass,Zeta0,Fre,mtmd,zetatmd,fretmd)

    omegatmd= 2*pi*fretmd;
    ktmd=mtmd.*omegatmd.^2;
    ctmd= 2.*mtmd.*omegatmd.*zetatmd;
    tmdnum=size(mtmd,1);
    disp("共设置"+num2str(tmdnum)+"个tmd");

    matrixsize=nModes+tmdnum;

      % Initialize the output matrix
      MM=zeros(matrixsize);
    
      for k1 = 1:matrixsize
          if k1==1
              MM(k1,k1)=Mass;
          elseif k1>=2
              MM(k1,k1)=mtmd(k1-1);
          end
      end
      clear k1
  
      CC  = zeros(size(MM, 1), size(MM, 2));
      CC1 = zeros(size(MM, 1), size(MM, 2));
      CC2 = zeros(size(MM, 1), size(MM, 2));
  
      for k1 = 1:matrixsize
  
              if k1<=1
                  CC1(k1,k1)=Zeta0 * 4 * pi * Mass * Fre;
              elseif k1>1
                  CC1(k1,k1)=ctmd(k1-1);
              end
      end
  
      for k1 = 1:matrixsize
          if k1<=1
              for k2 = 1:length(mtmd)
                  CC2(k1,k1)=CC2(k1,k1)+ctmd(k2);
              end
          elseif k1>1
              CC2(1,k1)=-ctmd(k1-1);
              CC2(k1,1)=-ctmd(k1-1);
          end
      end  
  
      CC = CC1 + CC2;
  
      clear k1
      KK = zeros(size(MM, 1), size(MM, 2));
      KK1 = zeros(size(MM, 1), size(MM, 2));
      KK2 = zeros(size(MM, 1), size(MM, 2));
      for k1 = 1:matrixsize
  
              if k1==1
                %   KK1(k1,k1)=Mass*(2*pi*Fre)^2-rho*U^2*H4;
                  KK1(k1,k1)=Mass*(2*pi*Fre)^2;
              elseif k1>1
                  KK1(k1,k1)=ktmd(k1-1);
              end
      end
      for k1 = 1:matrixsize
          if k1<=1
              for k2 = 1:length(mtmd)
                  KK2(k1,k1)=KK2(k1,k1)+ktmd(k2);
              end
          elseif k1>1
              KK2(1,k1)=-ktmd(k1-1);
              KK2(k1,1)=-ktmd(k1-1);
          end
      end
      KK = KK1 + KK2;
end
    