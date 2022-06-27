clc
clear
close all

%--------------------------------------------------------------------------
% Default inspection data
% Marra, A. M., et al. (2015). "Measurements and improved model of 
% vortex-induced vibration for an elongated rectangular cylinder." 
% Journal of Wind Engineering and Industrial Aerodynamics 147: 358-367.
% amplitude with Sc number
Sc_plot = [1.9;3.3;6.0;8.7;13.0;21.7;38.7;54.4;78.1];
Amplitude_plot = [272.8613;253.4274;211.9238;193.5358;161.001;123.0093;82.3407;58.9949;45.847]/291.8399*0.07;
PlotParameters = table(Sc_plot,Amplitude_plot);



%--------------------------------------------------------------------------
Fre_reference = [7.97;7.97;7.87;7.87;7.88;7.87;7.87;7.87;7.88];
Mass_reference = [6.99;6.99;7.19;7.19;7.19;7.19;7.19;7.19;7.19];
Zeta0_reference = [0.058;0.1;0.18;0.26;0.39;0.65;1.15;1.63;2.34]/100;
rho_reference = [1.19;1.19;1.22;1.22;1.22;1.22;1.21;1.22;1.22];
Sc_reference = [1.9;3.3;6.0;8.7;13.0;21.7;38.7;54.4;78.1];
DynamicParameters = table(Fre_reference,Mass_reference,Zeta0_reference,rho_reference,Sc_reference);



% aerodynamic parameters
Y1k_reference=[4.4557;5.0669;7.5116;8.2756;11.2552;17.5198;29.3616;41.0787;56.5112]/63.2728*60;
epsilonk_reference=[3.9198;4.0841;5.0597;4.8468;6.3544;9.5138;18.4515;35.5879;56.5112]/63.1228*12000;
Y2k_reference=[0;0;0;0;0;0;0;0;0];
ClK_reference=[0;0;0;0;0;0;0;0;0];
U_n0D_reference=[9.59501147074820;9.60967476438109;9.21413703087144;9.23010378861333;8.90889190943955;8.53309627387713;8.51366475455395;8.33081646873124;8.14314388632340];
AerodynamicParameters = table(Y1k_reference,epsilonk_reference,Y2k_reference,ClK_reference,U_n0D_reference);

collectdata=[];
for k1 =1:1
casenum=k1;
Fre=DynamicParameters.Fre_reference(casenum);
Mass=DynamicParameters.Mass_reference(casenum);
Zeta0=DynamicParameters.Zeta0_reference(casenum);
rho=DynamicParameters.rho_reference(casenum);
D=75/1000;
% D=100;
% D=0.5;
U_n0D=AerodynamicParameters.U_n0D_reference(casenum);
% U_n0D=mean(U_n0D_reference);
% U=Fre*D*U_n0D;
% U=20.92*4.185*55/1000;
U = 5.7354;



Y1k=AerodynamicParameters.Y1k_reference(casenum);
epsilonk=AerodynamicParameters.epsilonk_reference(casenum);
Y2k = AerodynamicParameters.Y2k_reference(casenum);
Clk = AerodynamicParameters.ClK_reference(casenum);
h=0.001;
T=80;
t=0:h:T;
P= zeros(1,length(t));
u0 = 0.05;  
udot0 = 0;

L=986/1000;
B=300/1000;
Sc=4*pi*Mass*Zeta0/rho/B/D/L;

out=NonlinearScanlan(Fre, Mass, Zeta0, rho, D, U, Y1k, epsilonk, Y2k, Clk, t, P, u0, udot0);

amplitude=max(out(2,round(end-end/5):end));

figure
% plot(out(1,round(end-end/30):end),out(2,round(end-end/30):end))
plot(out(1,1:end),out(2,1:end))

collectdata=[collectdata;[Sc amplitude]];
betacollectdata(k1)=sqrt(4/epsilonk*(1-Sc*2*pi*Fre*D/U*B/2/pi/Y1k/D));
end

% figure
% [psd_avg, f, psd_plot]=fft_transfer(1/h,out(2,1:end)');
% plot(f,psd_plot)
% 
% figure
% scatter(PlotParameters.Sc_plot,PlotParameters.Amplitude_plot,'blue','LineWidth',2)
% xlim([0 80])
% ylim([0 0.07])
% hold on
% scatter(collectdata(:,1),collectdata(:,2),'red','LineWidth',2)
% scatter(collectdata(:,1),betacollectdata,'green','LineWidth',2)
% legend('Experiment', 'Newmark beta method', 'Theory method');
% xlabel("$Sc$",'interpreter','latex')
% ylabel("$\sqrt(2)*y'/D$",'interpreter','latex')