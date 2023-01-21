%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: Shengyi xushengyichn@outlook.com
%Date: 2022-11-28 17:39:07
%LastEditors: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%LastEditTime: 2023-01-19 11:14:34
%FilePath: \NonlinearScanlan\20230116全阶模态对1TMD控制效果影响\data_analysis.m
%Description: 分析数据
%
%Copyright (c) 2022 by Shengyi xushengyichn@outlook.com, All Rights Reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 批量分析气动力作用于1-5阶模态时的临近模态影响

clc
clear 
close all

modes=1:5;
collectdata=[];
for k1 = modes
    % 读取数据
%     data1=importdata('10modes_onetmd_results_loc_mode1.mat');
    str1="data_mode=importdata('17modes_onetmd_results_loc_mode"+num2str(k1)+".mat');";
    eval(str1)
    str2="data_accurate=importdata('100modes_onetmd_results_loc_mode"+num2str(k1)+".mat');";
    eval(str2)
    clear str1 str2
    mode_index=data_mode.calmodes_index;
    nmodes_onetmd_results_loc=data_mode.nmodes_onetmd_results_loc;

    % 横坐标
    loc= data_accurate(:,1);


    % 提取考虑更多阶模态的位移
    str="dis_all_"+num2str(k1)+"=[];";
    eval(str)
    clear str
    
    for k2 = 1:length(mode_index)
        datasize=size(nmodes_onetmd_results_loc,1)/length(mode_index);
        str="dis_index"+num2str(k2)+"=nmodes_onetmd_results_loc(1+datasize*(k2-1):datasize*k2,3);";
        eval(str)
        clear str
        str="dis_all_"+num2str(k1)+"=[dis_all_"+num2str(k1)+" dis_index"+num2str(k2)+"];";
        eval(str)
        clear str
    end

    % 提取精确解位移
    dis_accurate=data_accurate(1:661,3);
   
    str="dis_accurate_"+num2str(k1)+"=data_accurate(1:661,3);";
    eval(str)
    clear str

    % 读取模态信息
    modeinfo = load('modeinfo_all.mat');
    nodegap=modeinfo.nodegap;
    mode=modeinfo.mode_re;
    for t01 = 1:length(loc)
        [~, index] = sort(abs(nodegap - loc(t01))); %查找与xTMD最接近的点的排序
        xResult = nodegap(index(1:2)); %获取最接近的两个点的x坐标
        mode2nodes = mode(index(1:2), 1:17); %获取两个点坐标的y值
        phi_result = interp1(xResult, mode2nodes, loc(t01), 'linear', 'extrap'); %插值以后任意点的振型
        mode_TMD_location(t01, 1:17) = phi_result(1:17);
    end
    
    
    % 计算相对精确解的位移差(每一阶的贡献值)
    for k2 = 1:length(mode_index)-1
%     for k2 = 2
        str="dis_"+num2str(mode_index(k2+1))+"_contri=(dis_index"+num2str(k2+1)+"-dis_index"+num2str(k2)+")./dis_accurate;";
        disp(str)
        eval(str)
        clear str
        mode_shape_control=mode_TMD_location(:,mode_index(1));%TMD所控制模态的振型
        mode_shape_consider=mode_TMD_location(:,mode_index(k2+1));
        % dis_1_plotdata=[dis_1_contri mode_shape_control mode_shape_consider dis_accurate]; %用于参考的代码，展示下方str的某一个特例
        % 第一列 某一个模态贡献值 第二列 TMD控制模态振型 第三列 考虑的模态振型 第四列 精确解位移
        % str="plotdata_dis_"+num2str(mode_index(k2+1))+"_for_mode_"+num2str(mode_index(1))+"=[dis_"+num2str(mode_index(k2+1))+"_contri mode_shape_control mode_shape_consider dis_accurate];";
        str="temp=[dis_"+num2str(mode_index(k2+1))+"_contri mode_shape_control mode_shape_consider dis_accurate];";
        eval(str)
        temp(:,1)=abs(temp(:,1));
        temp(:,2)=abs(temp(:,2))/max(abs(temp(:,2)));
        temp(:,3)=abs(temp(:,3))/max(abs(temp(:,3)));
        temp(:,4)=-normalize(temp(:,4),'range')+1;%缩放精确位移数据
        figure
        scatter(temp(:,2),temp(:,3),round(temp(:,4)*50,0)+10,temp(:,1))
%         scatter(temp(1:110,2),temp(1:110,3),round(temp(1:110,4)*50,0)+10,temp(1:110,1))
%         hold on
%         scatter(temp(330:440,2),temp(330:440,3),round(temp(330:440,4)*50,0)+10,temp(330:440,1))
        colorbar
        clear str
        str="plotdata_dis_"+num2str(mode_index(k2+1))+"_for_mode_"+num2str(mode_index(1))+"=temp;";
        eval(str)
        clear str
        str="collectdata.plotdata_dis_"+num2str(mode_index(k2+1))+"_for_mode_"+num2str(mode_index(1))+"=temp;";
        eval(str)
        clear str
    end
    collectdata.mode=mode_TMD_location;
    collectdata.loc=loc;
    str="collectdata.dis_all_"+num2str(k1)+"=dis_all_"+num2str(k1)+";";
    eval(str)
    clear str
    str="collectdata.dis_accurate_"+num2str(k1)+"=dis_accurate_"+num2str(k1)+";";
    eval(str)
    clear str
    save("mode_contribution_plotdata.mat","collectdata")
    % 画图代码
    % figure
    % scatter(abs(plotdata_dis_2_for_mode_1(:,2))/max(abs(plotdata_dis_2_for_mode_1(:,2))),abs(plotdata_dis_2_for_mode_1(:,3))/max(abs(plotdata_dis_2_for_mode_1(:,3))),plotdata_dis_2_for_mode_1(:,4)*100,plotdata_dis_2_for_mode_1(:,1))
    colorbar

%     figure
% %     
%     plot(loc,dis_accurate,'k','LineWidth',2)
%     hold on
% %     figure
%     for k2 = 1:length(mode_index)
%         str="plot(loc,dis_index"+num2str(k2)+");";
%         eval(str)
%         hold on
%     end
% % 
%     figure
%     for k2 = 1:length(mode_index)
%         plot(mode_TMD_location(:,k2))
%         hold on
%     end
% 
%         figure
%     for k2 = 1:3
%         str="plot(loc,dis_index"+num2str(k2)+");";
%         eval(str)
%         hold on
%     end
% 
%     figure
%     for k2 = 1:3
%         plot(mode_TMD_location(:,k2))
%         hold on
%     end
% 
%     figure
%     for k2 = 1:length(mode_index)-1
%         str="plot(loc,dis_"+num2str(mode_index(k2+1))+"_contri);";
%         eval(str)
%         hold on
%     end   
end


close all
return

%% 绘制分别考虑1-前5阶模态对TMD控制第一阶涡振效果的影响
clc; clear; close all;

% 读取数据
data = importdata('nmodes_onetmd_results_loc.mat');

loc=data(1:661,1);
dis=data(1:661,3);
figure
plot(loc,dis)

loc2=data(661*1+1:661*2,1);
dis2=data(661*1+1:661*2,3);
hold on
plot(loc2,dis2)

loc3=data(661*2+1:661*3,1);
dis3=data(661*2+1:661*3,3);
hold on
plot(loc3,dis3)

loc4=data(661*3+1:661*4,1);
dis4=data(661*3+1:661*4,3);
hold on
plot(loc4,dis4)

loc5=data(661*4+1:661*5,1);
dis5=data(661*4+1:661*5,3);
hold on
plot(loc5,dis5)

save nmodes_onetmd_results_loc_alldata.mat

%% 绘制考虑前8阶模态和考虑100阶的相对精确解的差异
clc
clear
close all

data = importdata('10modes_onetmd_results_loc.mat');
num=size(data,1)/9;
for k1 = 1:9
    str="dis"+num2str(data(1+num*(k1-1),2))+"=data(1+num*(k1-1):num*k1,3);";
    eval(str)
end
loc=data(1:num,1);

mode2effect=(dis2-dis1)./dis100;
mode3effect=(dis3-dis2)./dis100;
mode4effect=(dis4-dis3)./dis100;
mode5effect=(dis5-dis4)./dis100;
mode6effect=(dis6-dis5)./dis100;
mode7effect=(dis7-dis6)./dis100;
mode8effect=(dis8-dis7)./dis100;
mode100effect=(dis100-dis8)./dis100;

% figure
% plot(loc,mode2effect)
% hold on
% plot(loc,mode3effect)
% % plot(loc,mode4effect)
% % plot(loc,mode5effect)
% % plot(loc,mode6effect)
% % plot(loc,mode7effect)
% % plot(loc,mode8effect)
% % plot(loc,mode100effect)

% figure
% plot(loc,dis1)
% hold on
% plot(loc,dis2)
% plot(loc,dis3)
% % plot(loc,dis4)
% % plot(loc,dis5)
% % plot(loc,dis6)
% % plot(loc,dis7)
% % plot(loc,dis8)
% % plot(loc,dis100)

modeinfo = load('modeinfo_all.mat');
nodegap=modeinfo.nodegap;
mode=modeinfo.mode_re;
for t1 = 1:length(loc)

    [~, index] = sort(abs(nodegap - loc(t1))); %查找与xTMD最接近的点的排序
    xResult = nodegap(index(1:2)); %获取最接近的两个点的x坐标
    mode2nodes = mode(index(1:2), 1:8); %获取两个点坐标的y值
    phi_result = interp1(xResult, mode2nodes, loc(t1), 'linear', 'extrap'); %插值以后任意点的振型
    mode_TMD_location(t1, 1:8) = phi_result(1:8);
end
% figure
% plot(nodegap,mode(:,1))
% hold on
% plot(loc,modeTMD(:,1))

dis=[dis1 dis2 dis3 dis4 dis5 dis6 dis7 dis8 dis100];
modeeffect=[mode2effect mode3effect mode4effect mode5effect mode6effect mode7effect mode8effect];
Freq=modeinfo.Freq;

for k1 = 1:length(Freq)-1
    Freqeffect(k1,1)=(Freq(k1+1));
end



% 遍历所有点
figure
for k1 = 1:length(loc)
pointseq=find(loc==loc(k1));
point_mode_effect=modeeffect(pointseq,:);
point_freq_effect=Freqeffect(1:7);
point_mode_shape=mode_TMD_location(pointseq,:);

for k2 = 1:length(point_mode_shape)-1
%    point_mode_shape_effect(k1,1)=((abs(point_mode_shape(k1+1))-abs(point_mode_shape(1)))/point_mode_shape(1)*100);
   point_mode_shape_effect(k2,1)=(abs(point_mode_shape(k2+1))/(max(abs(mode_TMD_location(:,k2+1)))));
end
phi1=abs(point_mode_shape(1)/max(abs(mode_TMD_location(:,1)))*ones(7,1));
scatter3(point_mode_shape_effect,point_freq_effect,point_mode_effect,50*ones(7,1),phi1)
hold on
end
xlabel("mode shape")
ylabel("frequency")
zlabel("mode contribution")

colorbar


%% 按照模态来获取画图数据

clc
clear
close all

data = importdata('10modes_onetmd_results_loc.mat');
num=size(data,1)/9;
for k1 = 1:9
    str="dis"+num2str(data(1+num*(k1-1),2))+"=data(1+num*(k1-1):num*k1,3);";
    eval(str)
end
loc=data(1:num,1);

mode2contri=(dis2-dis1)./dis100;
mode3contri=(dis3-dis2)./dis100;
mode4contri=(dis4-dis3)./dis100;
mode5contri=(dis5-dis4)./dis100;
mode6contri=(dis6-dis5)./dis100;
mode7contri=(dis7-dis6)./dis100;
mode8contri=(dis8-dis7)./dis100;
mode100contri=(dis100-dis8)./dis100;

modecontri=[mode2contri mode3contri mode4contri mode5contri mode6contri mode7contri mode8contri];

for k1 = 1:7
    tempdata=zeros(size(loc,1),3);%获取绘图数据，第一列为模态1坐标，第二列为高阶模态坐标，第三列为高阶级模态贡献
    
    modeinfo = load('modeinfo_all.mat');
    nodegap=modeinfo.nodegap;
    mode=modeinfo.mode_re;
    for t01 = 1:length(loc)

        [~, index] = sort(abs(nodegap - loc(t01))); %查找与xTMD最接近的点的排序
        xResult = nodegap(index(1:2)); %获取最接近的两个点的x坐标
        mode2nodes = mode(index(1:2), 1:8); %获取两个点坐标的y值
        phi_result = interp1(xResult, mode2nodes, loc(t01), 'linear', 'extrap'); %插值以后任意点的振型
        mode_TMD_location(t01, 1:8) = phi_result(1:8);
    end

%     for t02 = 1:length(loc)
%         pointseq(t02)=find(loc==loc(t02));
%         point_mode_shape(t02)=modeTMD(pointseq,:);
%     end
    phi1=abs(mode_TMD_location(:,1)/max(abs(mode_TMD_location(:,1))));
    phin=abs(mode_TMD_location(:,k1+1)/max(abs(mode_TMD_location(:,k1+1))));
    tempdata(:,1)=phi1;
    tempdata(:,2)=phin;
    tempdata(:,3)=modecontri(:,k1);

    str="plotdata_mode"+num2str(k1+1)+"=tempdata;";
    eval(str)
end
figure
scatter(plotdata_mode2(:,1),plotdata_mode2(:,2),[],plotdata_mode2(:,3))
% scatter3(plotdata_mode2(:,1),plotdata_mode2(:,2),plotdata_mode2(:,3))
xlabel("modal displacement 1")
ylabel("modal displacement 2")
% zlabel("mode2 contribution value")

figure
scatter(plotdata_mode3(:,1),plotdata_mode3(:,2),[],plotdata_mode3(:,3))
% scatter3(plotdata_mode2(:,1),plotdata_mode2(:,2),plotdata_mode2(:,3))
xlabel("modal displacement 1")
ylabel("modal displacement 3")
% zlabel("mode2 contribution value")

figure
scatter(plotdata_mode5(:,1),plotdata_mode5(:,2),[],plotdata_mode5(:,3))
% scatter3(plotdata_mode2(:,1),plotdata_mode2(:,2),plotdata_mode2(:,3))
xlabel("modal displacement 1")
ylabel("modal displacement 5")
% zlabel("mode2 contribution value")


figure
scatter(plotdata_mode6(:,1),plotdata_mode6(:,2),[],plotdata_mode6(:,3))
% scatter3(plotdata_mode2(:,1),plotdata_mode2(:,2),plotdata_mode2(:,3))
xlabel("modal displacement 1")
ylabel("modal displacement 6")
% zlabel("mode2 contribution value")


%  f = fit([plotdata_mode2(:,1) plotdata_mode2(:,2)],plotdata_mode2(:,3),"poly11")
%  plot(f,[plotdata_mode2(:,1) plotdata_mode2(:,2)],plotdata_mode2(:,3))

