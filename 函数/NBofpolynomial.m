%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: Shengyi xushengyichn@outlook.com
%Date: 2022-05-15 16:16:48
%LastEditors: Shengyi xushengyichn@outlook.com
%LastEditTime: 2022-06-01 20:48:42
%FilePath: \NonlinearScanlan\NonlinearScanlan.m
%Description: 本函数目的为计算单自由度经验非线性模型响应，已经过验证计算无误。
%
%Copyright (c) 2022 by Shengyi xushengyichn@outlook.com, All Rights Reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a function
function out = NonlinearScanlan(Fre, Mass, Zeta0, rho, D, U, Y1k, epsilonk, Y2k, ClK, t, P, u0, udot0)

    % Nonlinear Newmark's Direct Integration Method with Nonlinear Scanlan
    % empirial noninear
    % model++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    % (n = number of time steps)
    % (ndof = number degrees of freedom)

    % INPUT
    % Fre      = Frequency of the system         =>[1]
    % Mass     = Mass of the system              =>[1]
    % Zeta0    = Damping ratio of the system     =>[1]
    % rho      = Air density                     =>[1]
    % D        = Reference length                =>[1]
    % U        = Wind speed at a certain reduced frequency       =>[1]
    % Y1k      = Y1 at a certain reduced frequency       =>[1]
    % epsilonk = Epsilon at a certain reduced frequency       =>[1]
    % Y2k      = Y2 at a certain reduced frequency       =>[1]
    % Clk      = Cl at a certain reduced frequency       =>[1]
    % t        = Time vector         =>[1,n]
    % P        = load vs. time       =>[1,n]
    % u0       = Initial displacements =>[1]
    % udot0    = Initial velocity =>[1]
    % gam      = gamma (constant)
    % beta     = beta  (constant)

    %--------------------------------------------------------------------------
    % beta = 0,     gamma = 1/2 => explicit central difference method
    % beta = 1/4,   gamma = 1/2 => undamped trapezoidal rule (implicit)

    %--------------------------------------------------------------------------
    % Default inspection data
    % Marra, A. M., et al. (2015). "Measurements and improved model of
    % vortex-induced vibration for an elongated rectangular cylinder."
    % Journal of Wind Engineering and Industrial Aerodynamics 147: 358-367.

    Fre_reference = [7.97; 7.97; 7.87; 7.87; 7.88; 7.87; 7.87; 7.87; 7.88];
    Mass_reference = [6.99; 6.99; 7.19; 7.19; 7.19; 7.19; 7.19; 7.19; 7.19];
    Zeta0_reference = [0.0058; 0.1; 0.18; 0.26; 0.39; 0.65; 1.15; 1.63; 2.34] / 100;
    rho_reference = [1.19; 1.19; 1.22; 1.22; 1.22; 1.22; 1.21; 1.22; 1.22];
    Sc_reference = [1.9; 3.3; 6.0; 8.7; 13.0; 21.7; 38.7; 54.4; 78.1];
    DynamicParameters = table(Fre_reference, Mass_reference, Zeta0_reference, rho_reference, Sc_reference);

    % aerodynamic parameters
    Y1k_reference = [4.4557; 5.0669; 7.5116; 8.2756; 11.2552; 17.5198; 29.3616; 41.0787; 56.5112] / 63.2728 * 60;
    epsilonk_reference = [4.0534; 4.2623; 5.0597; 4.8468; 6.3544; 9.5138; 18.4515; 35.5879; 56.5112] / 63.2728 * 12000;
    Y2k_reference = [0; 0; 0; 0; 0; 0; 0; 0; 0];
    ClK_reference = [0; 0; 0; 0; 0; 0; 0; 0; 0];
    AerodynamicParameters = table(Y1k_reference, epsilonk_reference, Y2k_reference, ClK_reference);

    % amplitude with Sc number
    Sc_plot = [1.9; 3.3; 6.0; 8.7; 13.0; 21.7; 38.7; 54.4; 78.1];
    Amplitude_plot = [96.8733; 89.5479; 76.0242; 69.2624; 57.6169; 43.7176; 29.4425; 20.2389; 14.9797] / 119.9763 * 0.08;
    PlotParameters = table(Sc_plot, Amplitude_plot);
    % scatter(PlotParameters.Sc_plot,PlotParameters.Amplitude_plot,'blue','LineWidth',2)
    % xlim([0 120])
    % ylim([0 0.08])

    gamma = 1/2; % Factor in the Newmark algorithm
    beta = 1/4; % Factor in the Newmark algorithm

    casenum = 1;
    %% Scanlan nonlinear model parameters
    if ~exist('Fre', 'var')
        Fre = DynamicParameters.Fre_reference(casenum);
        disp("Frequency of the system is using the default parameter: " + num2str(Fre))
    end

    if ~exist('Mass', 'var')
        Mass = DynamicParameters.Mass_reference(casenum);
        disp("Mass of the system is using the default parameter: " + num2str(Mass))
    end

    if ~exist('Zeta0', 'var')
        Zeta0 = DynamicParameters.Zeta0_reference(casenum);
        disp("Damping ratio of the system is using the default parameter: " + num2str(Zeta0))
    end

    if ~exist('rho', 'var')
        rho = DynamicParameters.rho_reference(casenum);
        disp("Air density is using the default parameter: " + num2str(rho))
    end

    if ~exist('D', 'var')
        % reference length
        D = 75/1000;
        disp("Reference length is using the default parameter: " + num2str(D))
    end

    if ~exist('U', 'var')
        %     U_n0D=9.5;
        % U=DynamicParameters.Fre(casenum)*D*U_n0D;
        U = 20.92 * 4.185 * 55/1000;
        disp("Wind speed at a certain reduced frequency is using the default parameter: " + num2str(U))
    end

    if ~exist('Y1k', 'var')
        Y1k = AerodynamicParameters.Y1k_reference(casenum);
        disp("Y1 at a certain reduced frequency is using the default parameter: " + num2str(Y1k))
    end

    if ~exist('epsilonk', 'var')
        epsilonk = AerodynamicParameters.epsilonk_reference(casenum);
        disp("Epsilon at a certain reduced frequency is using the default parameter: " + num2str(epsilonk))
    end

    if ~exist('Y2k', 'var')
        Y2k = AerodynamicParameters.Y2k_reference(casenum);
        disp("Y2 at a certain reduced frequency is using the default parameter: " + num2str(Y2k))
    end

    if ~exist('Clk', 'var')
        ClK = AerodynamicParameters.ClK_reference(casenum);
        disp("Cl at a certain reduced frequency is using the default parameter: " + num2str(ClK))
    end

    if ~exist('t', 'var')
        h = 0.001;
        T = 10;
        t = 0:h:T;
        disp("Calcute length of time is : " + num2str(T) + " seconds by default value. dt is" + num2str(h))
    else
        h = t(2) - t(1);
        T = t(end);
        disp("Calcute length of time is : " + num2str(T) + " seconds. dt is " + num2str(h))
    end

    if ~exist('P', 'var')
        P = zeros(1, length(t));
        pp = P;
        disp("Ignore the P by default.")
    else
        pp = P;
    end

    matrixsize = 1; %calculate one degree of freedom

    if ~exist('u0', 'var')
        u0 = zeros(matrixsize, 1);
        u0(1) = 0.05;
        disp("u0 is set to " + num2str(u0) + " by default")
    end

    if ~exist('udot0', 'var')
        udot0 = zeros(matrixsize, 1);
        udot0(1) = 0;
        disp("udot0 is set to " + num2str(udot0) + " by default")
    end

    Y1 = Y1k;
    epsilon = epsilonk;
    Y2 = Y2k;

    MM = Mass;
    CC = Zeta0 * 4 * pi * MM * Fre;
    KK = 4 * pi^2 * MM * Fre^2;

    %% Calculate the response

    nModes = 1;

    gfun = @(u, udot) Scanlan_nonlinear(u, udot, MM, CC, KK, gamma, beta, h, rho, U, D, Y1, epsilon, Y2, matrixsize, nModes);
    u = nonlinear_newmark_krenk(gfun, MM, pp, u0, udot0, gamma, beta, h);

    u_nounit = u / D;
    s = t * U / D;
%     out = [s; u_nounit];
    out = [t; u];
    %     figure
    %     plot(s,u_nounit)

    %     figure
    %     [psd_avg, f, psd_plot] = fft_transfer(1/h,u_nounit');
    %     plot(f,psd_plot)

    %% Nonlinear Newmark algorithm
    function [u, udot, u2dot] = nonlinear_newmark_krenk(gfun, MM, pp, u0, udot0, gamma, beta, h)
        % Initialize variables
        u = zeros(size(MM, 1), size(pp, 2));
        udot = zeros(size(MM, 1), size(pp, 2));
        u2dot = zeros(size(MM, 1), size(pp, 2));
        % Initial conditions
        u(:, 1) = u0; % Initial displacements
        udot(:, 1) = udot0; % Initial velocities
        g = gfun(u0, udot0); % Evaluate the nonlinear function
        u2dot(:, 1) = MM \ (pp(:, 1) - g); % Calculate the initial accelerations

        for ii = 1:size(pp, 2) - 1
            % Prediction step
            u2dot(:, ii + 1) = u2dot(:, ii); % Predicted accelerations
            udot(:, ii + 1) = udot(:, ii) + h * u2dot(:, ii); % Predicted velocities
            u(:, ii + 1) = u(:, ii) + h * udot(:, ii) + 1/2 * h^2 * u2dot(:, ii); % Predicted displacements
            Nit = 0; % Number of iterations
            konv = 0; % Has the iterations converged 0=No, 1=Yes
            % Reusidal calculation, system matrices and increment correction
            while Nit < 1000 && konv == 0
                Nit = Nit + 1;
                [g, Ks] = gfun(u(:, ii + 1), udot(:, ii + 1)); % Calculate function value and the tangent
                rr = pp(:, ii + 1) - MM * u2dot(:, ii + 1) - g; % Calculate residual
                du = Ks \ rr; % Increment correction
                u(:, ii + 1) = u(:, ii + 1) + du; % Add incremental correction to the displacement
                udot(:, ii + 1) = udot(:, ii + 1) + gamma * h / (beta * h^2) * du; % Incremental correction for velocities
                u2dot(:, ii + 1) = u2dot(:, ii + 1) + 1 / (beta * h^2) * du; % Incremental correction accelerations

                if sqrt(rr' * rr) / length(rr) < 1.0e-8 % Convergence criteria
                    konv = 1; % konv = 1 if the convergence criteria is fulfilled.
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
    % function [g, ks] = Scanlan_nonlinear(u, udot, MM, CC, KK, gamma, beta, h, rho, U, D, Y1, epsilon, Y2, matrixsize, mode, phiResultall, nodegap)
    %     diagMatrix = zeros(matrixsize, matrixsize);
    %     diagMatrix(mode, mode) = 1;
    %     % for k1 = 1:nModes
    %     %     diagMatrix(k1, k1) = 1;
    %     % end
        
    % if ~exist('phiResultall', 'var')
    %     phi2=1;
    %     phi4=1;
    %     disp("未输入振型，为单自由度问题。")
    % else
    %     phi=phiResultall(mode,:);

    %     for k1 =1:length(nodegap)-1
    %         dx(k1)=nodegap(k1+1)-nodegap(k1);
    %     end

    %     clear k1
    %     phi2=0;
    %     phi4=0;
    %     for k1 =1:length(dx)
    %         phi2=phi2+phi(k1)^2*dx(k1);
    %         phi4=phi4+phi(k1)^4*dx(k1);
    %     end
    %     clear k1
        
    % end


    %     g = (CC - diagMatrix * rho .* U .* D .* Y1 .* (phi2 - epsilon .* u.^2 .* phi4 ./ D.^2)) * udot + (KK) * u; % Function value
    %     kc = CC - diagMatrix * rho .* U .* D .* Y1 .* (phi2 - epsilon .* u.^2 .* phi4 ./ D.^2);
    %     u_udot = zeros(matrixsize, matrixsize);

    %     for k1 = 1:size(u_udot, 1)
    %         u_udot(k1, k1) = u(k1) * udot(k1);
    %     end

    %     kspring = KK + 2 .* rho .* U .* D .* Y1 .* epsilon .* phi4 ./ D^2 * u_udot;
    %     ks = kspring + gamma .* h ./ (beta .* h.^2) .* kc + 1 ./ (beta .* h.^2) .* MM; % Linearization
    % end

end
