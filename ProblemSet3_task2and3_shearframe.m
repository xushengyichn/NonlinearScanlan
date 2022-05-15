%% Initiate
clc
clear all
close all


%% Establish system
L=3.5;
I=180e-6;
E=210e9;
k=12*E*I/L^3
m=12e3;

M=diag([m m]);
K=k*[4 -2; -2 2];
C=3e-1*M+5e-4*K;
C(1,2)=C(1,2)+500;
%% Solve eigenvalue problem

[phi,lambda]=eig(K,M);
omega=lambda.^0.5;

psi=(phi.'*C*phi)./(2*omega)

%% Create time vector and load
T=100;
dt=0.01;
t=[0:dt:T];
nt=length(t);

X1=500*randn(1,nt);
X2=600*randn(1,nt);

figure(); hold on; grid on;
plot(t,X1,'b');
plot(t,X2,'r');
xlabel('Time [s]'); ylabel('Load [N]'); legend({'X1' 'X2'}); 

% plotscriptmain('h',8,'w',12,'name','Task_2b_1','titlesize',6,'labelsize',6,'ticksize',6,'path',cd,'open','yes');
xlim([2 5]);
% plotscriptmain('h',8,'w',12,'name','Task_2b_2','titlesize',6,'labelsize',6,'ticksize',6,'path',cd,'open','yes');
close all
%% Solve system response by time integration

u0=[0;0];
udot0=[0;0];
P=[X1;X2];
beta=0.25;
gam=0.5;

[y ydot y2dot] = NewmarkInt(t,M,C,K,P,gam,beta,u0,udot0);

figure(); hold on; grid on;
plot(t,y(1,:),'b');
plot(t,y(2,:),'r');
xlabel('Time [s]'); ylabel('Response [m]'); legend({'y1' 'y2'}); 

% plotscriptmain('h',8,'w',12,'name','Task_2c_1','titlesize',6,'labelsize',6,'ticksize',6,'path',cd,'open','yes');
xlim([20 30]);
% plotscriptmain('h',8,'w',12,'name','Task_2c_2','titlesize',6,'labelsize',6,'ticksize',6,'path',cd,'open','yes');
close all
%% Statistics of the response

std_y=std(y,0,2)
rho_y=corrcoef(y.')

%% Simulate multiple load realisations in a loop

nsim=50;

Y1=zeros(nsim,nt);
Y2=zeros(nsim,nt);

for j=1:nsim

X1=500*randn(1,nt);
X2=100*randn(1,nt);
P=[X1;X2];
[y ydot y2dot] = NewmarkInt(t,M,C,K,P,gam,beta,u0,udot0);

%Collect the response of y1 and y2 the matrices Y1 and Y2
Y1(j,:)=y(1,:);
Y2(j,:)=y(2,:);

end

% Fourier transform of each of the response time series from the
% simulations
[f,G_Y1]=fft_function(Y1,dt);
[f,G_Y2]=fft_function(Y2,dt);

% Estimate the spectra directly from time series using the DFT
S_Y1=G_Y1.*conj(G_Y1)*T;
S_Y2=G_Y2.*conj(G_Y2)*T;

%Take the average of all specta
S_y1y1=mean(S_Y1,1);
S_y2y2=mean(S_Y2,1);

%Plot the spectra
figure(); hold on; grid on;
plot(f,S_y1y1,'b');
plot(f,S_y2y2,'r');
xlabel('Frequency [Hz]'); ylabel('Spectral density [m^2/Hz]'); legend({'S_{y1y1}' 'S_{y2y2}'}); 


xlim([-12 12]);
% plotscriptmain('h',8,'w',12,'name','Task_3a','titlesize',6,'labelsize',6,'ticksize',6,'path',cd,'open','yes');
%plotscriptmain('h',8,'w',12,'name','Task_3b','titlesize',6,'labelsize',6,'ticksize',6,'path',cd,'open','yes');
close all

%% Plot transfer function

for k=1:length(f)
    om=f(k)*2*pi;
    H(:,:,k)=inv(-om.^2*M+i*om*C+K);
    H_abs(:,:,k)= H(:,:,k)*H(:,:,k)';

end

figure(); hold on; grid on;
plot(f*2*pi,squeeze(H_abs(1,1,:)),'b');
plot(f*2*pi,squeeze(H_abs(2,2,:)),'r');

