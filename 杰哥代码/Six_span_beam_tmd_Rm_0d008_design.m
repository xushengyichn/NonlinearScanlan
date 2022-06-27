%%%VIV control of a 6-span beam
clc;clear all;close all
tic
%%%mode shape of a 6-span beam
lamda=[3.141592654	3.260534994	3.55640846	3.926602312	4.297529693	4.601417878];
z_virtue=[0.01:0.01:1];z_virtue=z_virtue';
z=[0:0.01:6];z=z';

for ii=1:6
    A(1)=1;B(1)=0;
    mode_virtue(:,1)=A(1)*(sinh(lamda(ii))*sin(lamda(ii)*z_virtue)-sin(lamda(ii))*sinh(lamda(ii)*z_virtue))+B(1)*(sinh(lamda(ii))*sin(lamda(ii)*(1-z_virtue))-sin(lamda(ii))*sinh(lamda(ii)*(1-z_virtue)));
    
    for jj=2:6
        A(jj)=2*A(jj-1)*(sinh(lamda(ii))*cos(lamda(ii))-sin(lamda(ii))*cosh(lamda(ii)))/(sinh(lamda(ii))-sin(lamda(ii)))-B(jj-1);
        B(jj)=A(jj-1);
        mode_virtue(:,jj)=A(jj)*(sinh(lamda(ii))*sin(lamda(ii)*z_virtue)-sin(lamda(ii))*sinh(lamda(ii)*z_virtue))+B(jj)*(sinh(lamda(ii))*sin(lamda(ii)*(1-z_virtue))-sin(lamda(ii))*sinh(lamda(ii)*(1-z_virtue)));
    end
    mode(:,ii)=[mode_virtue(:,1);mode_virtue(:,2);mode_virtue(:,3);mode_virtue(:,4);mode_virtue(:,5);mode_virtue(:,6);];
    mode(:,ii)=mode(:,ii)/max(mode(:,ii));
end

%%%modal mass and frequency
load prototype
mode(:,7)=ones(1,600);
m=26e3;dl=900/600;
for ii=1:7
    m_modal(ii)=sum(m*mode(:,ii).^2)*dl;
    integral_2(ii)=sum(mode(:,ii).^2)*dl;
    integral_3(ii)=sum(abs(mode(:,ii)).^3)*dl;
    integral_4(ii)=sum(mode(:,ii).^4)*dl;
    integral_5(ii)=sum(abs(mode(:,ii)).^5)*dl;
    integral_6(ii)=sum(mode(:,ii).^6)*dl;
    integral_7(ii)=sum(abs(mode(:,ii)).^7)*dl;
end
frequency_modal = [0.43 0.47 0.55 0.68 0.81 0.93 1];

mode_number=1;
rou=1.225;%%主梁频率、阻尼、质量、空气密度、风速
AERO_DATA=[8.144 -132 34516.7273160516 -2927803.46063482 108484017.986453 -1803274448.49376 10992199274.0953];
Ur=AERO_DATA(1);K1=2*pi/Ur;a1=AERO_DATA(2);a2=AERO_DATA(3);a3=AERO_DATA(4);a4=AERO_DATA(5);a5=AERO_DATA(6);a6=AERO_DATA(7);

L=150*6;l=1/dl:dl:L;%%主梁长度
eps=0.0/100;omg=frequency_modal(mode_number)*2*pi;B=20;D=B/4;
%phi=ones(length(l),1);
phi=mode(:,mode_number);
mode_integral_2=integral_2(mode_number);
mode_integral_3=integral_3(mode_number);
mode_integral_4=integral_4(mode_number);
mode_integral_5=integral_5(mode_number);
mode_integral_6=integral_6(mode_number);
mode_integral_7=integral_7(mode_number);
m_bridge=m_modal(mode_number);
c_bridge=m_bridge*2*eps*omg;k_bridge=m_bridge*omg^2;

k_bridge=k_bridge;c_bridge=c_bridge;
U=Ur*sqrt(k_bridge/m_bridge)/2/pi*D;

%NES postition
NES_position=50;%%%change of damper position
mode_factor=phi(NES_position,:);

%%%TMD parameters
m_damper=0.008*m*L;
for iii=201:201
    omg_damper=(0.8+(iii-1)*0.002)*omg;
    k_damper=m_damper*(omg_damper)^2;
    for jjj=201:201
        eps_damper=0.02+0.001*(jjj-1);
        c_damper=m_damper*2*eps_damper*omg_damper;
        
        MM = diag([m_bridge m_damper]); % Mass matrix
        h = 0.1; % Time step
        t = 0:h:500; % Time
        p = zeros(2,length(t)); % Initialize external load
        %p(1,:)=[timeseries' 0];
        gamma = 1/2; % Parameter in the Newmark algorithm
        beta = 1/4; % Parameter in the Newmark algorithm
        gfun = @(u,udot) bridge_damper(u,udot,rou,U,D,a1,a2,a3,a4,a5,a6,c_bridge,k_bridge,c_damper,k_damper,mode_factor,mode_integral_2,mode_integral_3,mode_integral_4,mode_integral_5,mode_integral_6,mode_integral_7,MM,gamma,beta,h); % Handle to the nonlinear function
        u0 = [0.01 0]'; % Initial displacement CAREFUL
        udot0 = [0 0]'; % Initial velocity
        [u, udot, u2dot] = nonlinear_newmark_krenk(gfun,MM,p,u0,udot0,gamma,beta,h); % Solve the response by the Nonlinear Newmark algorithm
        % Plot the response
        u=u';u2dot=u2dot';
        
        amp_dis_bridge(iii,jjj)=max(u((length(u)-1)*8/10:length(u),1));
        amp_dis_damper(iii,jjj)=max(u((length(u)-1)*8/10:length(u),2));
        std_dis_bridge(iii,jjj)=std(u((length(u)-1)*8/10:length(u),1));
        std_dis_damper(iii,jjj)=std(u((length(u)-1)*8/10:length(u),2));
        amp_acc_bridge(iii,jjj)=max(u2dot((length(u2dot)-1)*8/10:length(u2dot),1));
        amp_acc_damper(iii,jjj)=max(u2dot((length(u2dot)-1)*8/10:length(u2dot),2));
        std_acc_bridge(iii,jjj)=std(u2dot((length(u2dot)-1)*8/10:length(u2dot),1));
        std_acc_damper(iii,jjj)=std(u2dot((length(u2dot)-1)*8/10:length(u2dot),2));
        figure;plot(u(:,1));%hold on;plot(u(:,2),'r')
    end
end
toc

save tmd_0d008_amp_dis_bridge amp_dis_bridge
save tmd_0d008_amp_dis_damper amp_dis_damper
save tmd_0d008_std_dis_bridge std_dis_bridge
save tmd_0d008_std_dis_damper std_dis_damper
save tmd_0d008_amp_acc_bridge amp_acc_bridge
save tmd_0d008_amp_acc_damper amp_acc_damper
save tmd_0d008_std_acc_bridge std_acc_bridge
save tmd_0d008_std_acc_damper std_acc_damper
min_response=[min(min(amp_dis_bridge)) min(min(amp_dis_damper)) min(min(std_dis_bridge)) min(min(std_dis_damper)) min(min(amp_acc_bridge)) min(min(amp_acc_damper)) min(min(std_acc_bridge)) min(min(std_acc_damper))];

%Result=[amp_bridge(11,11) amp_damper(11,11)]
% [ envelope_bridge,fre ] = ee(u(:,1),h );
% [ envelope_damper,fre ] = ee(u(:,2),h );

% figure('color',[1 1 1]);
% %contour(m_damper*2*0.01*omg:m_damper*2*0.02*omg:m_damper*2*0.59*omg,0.1*m_bridge*omg^2:0.1*m_bridge*omg^2:4*m_bridge*omg^2,amp_bridge,5,'ShowText','on')
% [MMM,ccc]= contour(m_damper*2*0.01*omg:m_damper*2*0.02*omg:m_damper*2*0.59*omg,0.1*m_bridge*omg^2:0.1*m_bridge*omg^2:4*m_bridge*omg^2,amp_bridge,'ShowText','on','LevelList',[0.04 0.06 0.08 0.10]);
% ccc.LineWidth = 3;
% xlim([0.5e5 3.5e5]);ylim([1e8 7e8])
% xlabel('{\itc}_{\itd}');ylabel('{\itk}_{\itd}');%set(gca,'YTick',0.8:0.1:1.2);
% set(gca,'LineWidth',1.3,'FontName','Times New Roman','FontSize',15);
% clabel(MMM,ccc,'FontSize',12,'FontName','Times New Roman')
% axes('Position',[0.2 0.6 0.2 0.2]);
% contour(m_damper*2*0.01*omg:m_damper*2*0.02*omg:m_damper*2*0.59*omg,0.1*m_bridge*omg^2:0.1*m_bridge*omg^2:4*m_bridge*omg^2,amp_bridge,'ShowText','on','LevelList',[0.04 0.06 0.08 0.10]);
% xlim([2e5 2.5e5]);ylim([5e8 6e8])

%% Function file for bridge-damper system
function [g,ks]=bridge_damper(u,udot,rou,U,D,a1,a2,a3,a4,a5,a6,c_bridge,k_bridge,c_damper,k_damper,mode_factor,mode_integral_2,mode_integral_3,mode_integral_4,mode_integral_5,mode_integral_6,mode_integral_7,MM,gamma,beta,h)
g1=c_bridge*udot(1,1)+k_bridge*u(1,1)+c_damper*(udot(1,1)*mode_factor-udot(2,1))*mode_factor+k_damper*(u(1,1)*mode_factor-u(2,1))*mode_factor+rou*U^2*D*(a1*mode_integral_2+a2*abs(u(1,1)/D)*mode_integral_3+a3*(u(1,1)/D)^2*mode_integral_4+a4*abs(u(1,1)/D)^3*mode_integral_5+a5*(u(1,1)/D)^4*mode_integral_6+a6*abs(u(1,1)/D)^5*mode_integral_7)*udot(1,1)/U;
g2=c_damper*(udot(2,1)-udot(1,1)*mode_factor)+k_damper*(u(2,1)-u(1,1)*mode_factor);
g = [g1 g2]'; % Nonlinear function value
ct=[c_bridge+c_damper*mode_factor*mode_factor+rou*U^2*D*(a1/U*mode_integral_2+a2/U*abs(u(1,1)/D)*mode_integral_3+a3/U*(u(1,1)/D)^2*mode_integral_4+a4/U*abs(u(1,1)/D)^3*mode_integral_5+a5/U*(u(1,1)/D)^4*mode_integral_6+a6/U*abs(u(1,1)/D)^5*mode_integral_7) -c_damper*mode_factor; -c_damper*mode_factor c_damper];
kt=[k_bridge+k_damper*mode_factor*mode_factor+rou*U^2*D*(a2*sign(u(1,1))/D*mode_integral_3+a3*2*(u(1,1)/D)/D*mode_integral_4+a4*3*sign(u(1,1))*(u(1,1)/D)^2/D*mode_integral_5+a5*4*(u(1,1)/D)^3/D*mode_integral_6+a6*5*sign(u(1,1))*(u(1,1)/D)^4/D*mode_integral_7)*udot(1,1)/U -k_damper*mode_factor; -k_damper*mode_factor k_damper];
ks = kt +gamma/(beta*h)*ct+MM*1/(beta*h^2); % Linearization
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