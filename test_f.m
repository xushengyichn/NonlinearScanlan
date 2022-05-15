clc
clear
close all
% This is a function 
% function out=test_f(arg1, arg2, arg3)


% Nonlinear Newmark's Direct Integration Method with Nonlinear Scanlan
% empirial noninear model
%--------------------------------------------------------------------------
% OUTPUT
% u =       Displacemente Response   [n,1]
% (n = number of time steps)
% (ndof = number degrees of freedom)

% INPUT
% t =       Time vector         [1,n]
% M =       mass matrix         [ndof,ndof]
% C =       damping matrix      [ndof,ndof]
% K =       stiffness matrix    [ndof,ndof]
% P =       load vs. time       [ndof,n]
% gam =     gamma (constant)
% beta =    beta  (constant)
%u0 =       Initial displacements
%udot0 =    Initial velocity

%--------------------------------------------------------------------------
% beta = 0,     gamma = 1/2 -> explicit central difference method
% beta = 1/4,   gamma = 1/2 -> undamped trapezoidal rule (implicit)

%--------------------------------------------------------------------------
% Default inspection data
% Marra, A. M., et al. (2015). "Measurements and improved model of 
% vortex-induced vibration for an elongated rectangular cylinder." 
% Journal of Wind Engineering and Industrial Aerodynamics 147: 358-367.

Fre = [7.97;7.97;7.87;7.87;7.88;7.87;7.87;7.87;7.88];
Mass = [6.99;6.99;7.19;7.19;7.19;7.19;7.19;7.19;7.19];
Zeta0 = [0.0058;0.1;0.18;0.26;0.39;0.65;1.15;1.63;2.34]/100;
rho = [1.19;1.19;1.22;1.22;1.22;1.22;1.21;1.22;1.22];
Sc = [1.9;3.3;6.0;8.7;13.0;21.7;38.7;54.4;78.1];

DynamicParameters = table(Fre,Mass,Zeta0,rho,Sc);

% reference length
D = 75/1000; 

% aerodynamic parameters
Y1k=[4.4557;5.0669;7.5116;8.2756;11.2552;17.5198;29.3616;41.0787;56.5112]/63.2728*60;
epsilonk=[4.0534;4.2623;5.0597;4.8468;6.3544;9.5138;18.4515;35.5879;56.5112]/63.2728*12000;
Y2k=[0;0;0;0;0;0;0;0;0];
ClK=[0;0;0;0;0;0;0;0;0];
AerodynamicParameters = table(Y1k,epsilonk,Y2k,ClK);

% amplitude with Sc number
Sc_plot = [1.9;3.3;6.0;8.7;13.0;21.7;38.7;54.4;78.1];
Amplitude_plot = [96.8733;89.5479;76.0242;69.2624;57.6169;43.7176;29.4425;20.2389;14.9797]/119.9763*0.08;
PlotParameters = table(Sc_plot,Amplitude_plot);
% scatter(PlotParameters.Sc_plot,PlotParameters.Amplitude_plot,'blue','LineWidth',2)
% xlim([0 120])
% ylim([0 0.08])

gamma = 1/2; % Factor in the Newmark algorithm 
beta = 1/4; % Factor in the Newmark algorithm 

%     if ~exist('arg1', 'var')
%         arg1 = 1;
%         disp("arg1 is using the default parameter: "+num2str(arg1))
%     end
%     if ~exist('arg2', 'var')
%         arg2 = 2;
%         disp("arg2 is using the default parameter: "+num2str(arg1))
%     end
%     if ~exist('arg3', 'var')
%         arg3 = 3;
%         disp("arg3 is using the default parameter: "+num2str(arg1))
%     end
%     out=arg1+arg2+arg3;
    datanum=6;
    U_n0D=9.5;
    U=DynamicParameters.Fre(datanum)*D*U_n0D;
    %% Scanlan nonlinear model parameters
    rho=DynamicParameters.rho(datanum);
    U=20.92*4.185*55/1000;
    D=D;
    Y1=AerodynamicParameters.Y1k(datanum);
    epsilon=AerodynamicParameters.epsilonk(datanum);
    Y2=AerodynamicParameters.Y2k(datanum);

    MM=DynamicParameters.Mass(datanum);
    CC=DynamicParameters.Zeta0(datanum)*4*pi*MM*DynamicParameters.Fre(datanum);
    KK=4*pi^2*MM*DynamicParameters.Fre(datanum)^2;
    
    
    
    
    %% Calculate the response
    h=0.001;
    T=10;
    t=0:h:T;
    P0=0;

    matrixsize=1;
    nModes=1;
    nbeta=0.25;
    ngam=0.5;
    u0=zeros(matrixsize,1);
    u0(1)=0.05;
    udot0=zeros(matrixsize,1);
    udot0(1)=0;
    pp=zeros(1,length(t));
    gfun = @(u,udot) Scanlan_nonlinear(u,udot,MM,CC,KK,gamma,beta,h,rho,U,D,Y1,epsilon,Y2,matrixsize,nModes);
    u = nonlinear_newmark_krenk(gfun,MM,pp,u0,udot0,gamma,beta,h);
    
    u_nounit=u/D;
    s=t*U/D;
    figure
    plot(s,u_nounit)

%     figure
%     [psd_avg, f, psd_plot] = fft_transfer(1/h,u_nounit');
%     plot(f,psd_plot)
    
    
    
    
    
    
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
    
        g = (CC-diagMatrix*rho.*U.*D.*Y1.*(1-epsilon.*u.^2./D.^2))*udot+(KK)*u; % Function value
        kc= CC-diagMatrix*rho.*U.*D.*Y1.*(1-epsilon.*u.^2./D.^2);
        u_udot=zeros(matrixsize,matrixsize);
        for k1=1:size(u_udot,1)
            u_udot(k1,k1)=u(k1)*udot(k1);
        end
        kspring= KK+2.*rho.*U.*D.*Y1.*epsilon./D^2*u_udot;
        ks = kspring + gamma.*h./(beta.*h.^2).*kc + 1./(beta.*h.^2).*MM; % Linearization
    end
% end
