%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: Shengyi Xu xushengyichn@outlook.com
%Date: 2022-07-07 20:50:55
%LastEditors: Shengyi Xu xushengyichn@outlook.com
%LastEditTime: 2022-07-07 20:57:20
%FilePath: \twindeck_ID\polynomial_NB_adstiff.m
%Description: 本函数目的为计算多项式模型的响应，同时考虑气动刚度H4
%
%Copyright (c) 2022 by Shengyi Xu xushengyichn@outlook.com, All Rights Reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a function
function out = polynomial_NB_adstiff(Fre, Mass, Zeta0, rho, D, U, a,H4, t, P, u0, udot0)

    % Nonlinear Newmark's Direct Integration Method with polynomial model
    % (n = number of time steps)
    % (ndof = number degrees of freedom)

    % INPUT
    % Fre      = Frequency of the system         =>[1]
    % Mass     = Mass of the system              =>[1]
    % Zeta0    = Damping ratio of the system     =>[1]
    % rho      = Air density                     =>[1]
    % D        = Reference length                =>[1]
    % U        = Wind speed at a certain reduced frequency       =>[1]
    % a        = 多项式模型系数（四次多项式5个系数）       =>[1,5]
    % H4       = 气动刚度项目       =>[1,1]
    % t        = Time vector         =>[1,n]
    % P        = load vs. time       =>[1,n]
    % u0       = Initial displacements =>[1]
    % udot0    = Initial velocity =>[1]
    % gam      = gamma (constant)
    % beta     = beta  (constant)

    %--------------------------------------------------------------------------
    % beta = 0,     gamma = 1/2 => explicit central difference method
    % beta = 1/4,   gamma = 1/2 => undamped trapezoidal rule (implicit)

    h = t(2) - t(1);

    gamma = 1/2; % Factor in the Newmark algorithm
    beta = 1/4; % Factor in the Newmark algorithm

    a1 = a(1);
    a2 = a(2);
    a3 = a(3);
    a4 = a(4);
    a5 = a(5);
    %     a=zeros(5,1);
    %     a1=a(1);
    %     a2=a(2);
    %     a3=a(3);
    %     a4=a(4);
    %     a5=a(5);

    m = Mass;
    omega0 = 2 * pi * Fre;
    %     b1=rho*U*D*a1/m-2*Zeta0*omega0;
    b1 = rho * U * D * a1 / m;
    b2 = rho * U * a2 / m;
    b3 = rho * U * a3 / m / D;
    b4 = rho * U * a4 / m / D^2;
    b5 = rho * U * a5 / m / D^3;

    b = [b1 b2 b3 b4 b5] .* m;
    
    matrixsize = 1; %calculate one degree of freedom
    
    

    MM = Mass;
    CC = Zeta0 * 4 * pi * MM * Fre;
    % KK = 4 * pi^2 * MM * Fre^2;
    KK =MM*((2*pi*Fre)^2-rho*U^2*H4/MM);
    pp = P;
    %% Calculate the response
    y = 0:0.001:0.1;
    z = -b(1) - b(2) .* abs(y) - b(3) .* y.^2 - b(4) * abs(y).^3 - b(5) .* y.^4;
    figure
    plot(y, z)

    nModes = 1;

    gfun = @(u, udot) polynomial_nonlinear(u, udot, MM, CC, KK, gamma, beta, h, b, matrixsize, nModes);
    u = nonlinear_newmark_krenk(gfun, MM, pp, u0, udot0, gamma, beta, h)';

    % u_nounit = u / D;
    % s = t * U / D;
    %     out = [s; u_nounit];
    out = [t u];
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

                if sqrt(rr' * rr) / length(rr) < 1.0e-4 % Convergence criteria
                    konv = 1; % konv = 1 if the convergence criteria is fulfilled.
                end

                if Nit == 999 % Convergence criteria
                    error("迭代未收敛，此时为第" + num2str(ii) + "个循环")
                end

            end

        end

    end

    %% Function file for the Scanlan nonlinear model
    %     function [g,ks]= polynomial_nonlinear(u, udot, MM, CC, KK, gamma, beta, h, b, matrixsize, nModes)
    %         diagMatrix=zeros(matrixsize,matrixsize);
    %         for k1=1:nModes
    %             diagMatrix(k1,k1)=1;
    %         end
    %
    %         g = (CC-diagMatrix.*b(5).*u.^4-diagMatrix.*b(4).*abs(u).^3-diagMatrix.*b(3).*u.^2-diagMatrix.*b(2).*abs(u)-b(1))*udot+(KK).*u; % Function value
    %         kc= CC-diagMatrix.*b(5).*u.^4-diagMatrix.*b(4).*abs(u).^3-diagMatrix.*b(3).*u.^2-diagMatrix.*b(2).*abs(u)-b(1);
    %         kspring= KK+diagMatrix.*udot.*(-4.*b(5).*u^3-3*b(4)*u^3./abs(u)-2.*b(3).*u-b(2).*u./abs(u));
    %         ks = kspring + gamma.*h./(beta.*h.^2).*kc + 1./(beta.*h.^2).*MM; % Linearization
    %     end

    %% Function file for the polynomial model
    function [g, ks] = polynomial_nonlinear(u, udot, MM, CC, KK, gamma, beta, h, b, matrixsize, nModes)

        g = (CC - b(5) .* u.^4 - b(4) .* abs(u).^3 - b(3) .* u.^2 - b(2) .* abs(u) - b(1)) * udot + (KK) .* u; % Function value
        kc = CC - b(5) .* u.^4 - b(4) .* abs(u).^3 - b(3) .* u.^2 - b(2) .* abs(u) - b(1);
        kspring = KK + udot .* (-4 .* b(5) .* u^3 - 3 * b(4) * u^3 ./ abs(u) - 2 .* b(3) .* u - b(2) .* u ./ abs(u));
        ks = kspring + gamma .* h ./ (beta .* h.^2) .* kc + 1 ./ (beta .* h.^2) .* MM; % Linearization
    end

end
