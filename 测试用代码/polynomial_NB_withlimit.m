%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: Shengyi xushengyichn@outlook.com
%Date: 2022-05-15 16:16:48
%LastEditors: Shengyi Xu xushengyichn@outlook.com
%LastEditTime: 2022-07-11 10:12:41
%FilePath: \twindeck_ID\polynomial_NB_withlimit.m
%Description:
%本函数目的为计算多项式模型的响应。采用上下限设置，使得不同振幅使用不同的函数，但是，由于存在斜率突变点，非常容易不收敛，计算结果可能也是有误的，因此决定放弃。2022年7月11日又发现可能是由于振幅的计算错误导致，原先使用的是位移，实际应该是sqrt(u^2+(udot/omega)^2),需要重新验证
%
%Copyright (c) 2022 by Shengyi xushengyichn@outlook.com, All Rights Reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a function
function out = polynomial_NB_withlimit(Fre, Mass, Zeta0, rho, D, U, a, upperlimit, lowerlimit, t, P, u0, udot0)

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
    % upperlimit        = 所能涵盖振幅最大值       =>[1]
    % lowerlimit        = 所能涵盖振幅最小值       =>[1]
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

    a1_lower=a1+4/3*a2/pi*lowerlimit+a3/4*lowerlimit^2+8/15*a4/pi*lowerlimit^3+a5/8*lowerlimit^4;
    a1_upper=a1+4/3*a2/pi*upperlimit+a3/4*upperlimit^2+8/15*a4/pi*upperlimit^3+a5/8*upperlimit^4;
    a1_lower=3.7358;
    a1_upper=1.9081e+41;
    m = Mass;
    omega0 = 2 * pi * Fre;
    %     b1=rho*U*D*a1/m-2*Zeta0*omega0;
    b1 = rho * U * D * a1 / m;
    b2 = rho * U * a2 / m;
    b3 = rho * U * a3 / m / D;
    b4 = rho * U * a4 / m / D^2;
    b5 = rho * U * a5 / m / D^3;

    b1_lower=rho*U*D*a1_lower/m;
    b1_upper=rho*U*D*a1_upper/m;
    b_lower = [b1_lower 0 0 0 0] .*m;
    b_upper = [b1_upper 0 0 0 0] .*m;

    b = [b1 b2 b3 b4 b5] .* m;

    lowerlimit=lowerlimit*D;%Newmark beta法中采用实际位移
    upperlimit=upperlimit*D;

    
    

%     b_lower=b1+b2*abs(lowerlimit)+b3*lowerlimit^2+b4*abs(lowerlimit)^3+b5*lowerlimit^4;%直接用启动阻尼等效来计算b的参数似乎有问题
% 
%     b_upper=b1+b2*abs(upperlimit)+b3*upperlimit^2+b4*abs(upperlimit)^3+b5*upperlimit^4;

    matrixsize = 1; %calculate one degree of freedom

    MM = Mass;
    CC = Zeta0 * 4 * pi * MM * Fre;
    KK = 4 * pi^2 * MM * Fre^2;
    pp = P;
    %% Calculate the response
%     y = 0:0.001:0.1;
%     z = -b(1) - b(2) .* abs(y) - b(3) .* y.^2 - b(4) * abs(y).^3 - b(5) .* y.^4;
%     figure
%     plot(y, z)

    nModes = 1;

    gfun1 = @(u, udot) polynomial_nonlinear(u, udot, MM, CC, KK, gamma, beta, h, b, matrixsize, nModes);
    gfun2 = @(u, udot) polynomial_nonlinear(u, udot, MM, CC, KK, gamma, beta, h, b_lower, matrixsize, nModes);
    gfun3 = @(u, udot) polynomial_nonlinear(u, udot, MM, CC, KK, gamma, beta, h, b_upper, matrixsize, nModes);
    % gfun2 = @(u, udot) limitdampingratio(u, udot, MM, CC, KK, gamma, beta, h, lowerlimit, matrixsize, nModes);
    % gfun3 = @(u, udot) limitdampingratio(u, udot, MM, CC, KK, gamma, beta, h, upperlimit, matrixsize, nModes);

    u = nonlinear_newmark_krenk(gfun1,gfun2,gfun3, MM, pp, u0, udot0, gamma, beta, h, upperlimit,lowerlimit)';

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
    function [u, udot, u2dot] = nonlinear_newmark_krenk(gfun1,gfun2,gfun3, MM, pp, u0, udot0, gamma, beta, h,upperlimit,lowerlimit)
        % Initialize variables
        u = zeros(size(MM, 1), size(pp, 2));
        udot = zeros(size(MM, 1), size(pp, 2));
        u2dot = zeros(size(MM, 1), size(pp, 2));
        % Initial conditions
        u(:, 1) = u0; % Initial displacements
        udot(:, 1) = udot0; % Initial velocities
        if abs(u(:, 1)) < lowerlimit
            g = gfun2(u(:, 1), udot0);
        else
            if abs(u(:, 1)) > upperlimit
                g = gfun3(u(:, 1), udot0);
                disp('gfun3')
            else
                g = gfun1(u(:, 1), udot0);
                disp('gfun1')
            end
        end
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
                if abs(u(:, ii + 1)) < lowerlimit
                    [g, Ks] = gfun2(u(:, ii + 1), udot(:, ii + 1));
                else
                    if abs(u(:, ii + 1)) > upperlimit
                        [g, Ks] = gfun3(u(:, ii + 1), udot(:, ii + 1));
                        disp('gfun3')
                    else
                        [g, Ks] = gfun1(u(:, ii + 1), udot(:, ii + 1));
                        disp('gfun1')
                    end
                end
                % [g, Ks] = gfun(u(:, ii + 1), udot(:, ii + 1)); % Calculate function value and the tangent
                rr = pp(:, ii + 1) - MM * u2dot(:, ii + 1) - g; % Calculate residual
                du = Ks \ rr; % Increment correction
                u(:, ii + 1) = u(:, ii + 1) + du; % Add incremental correction to the displacement
                udot(:, ii + 1) = udot(:, ii + 1) + gamma * h / (beta * h^2) * du; % Incremental correction for velocities
                u2dot(:, ii + 1) = u2dot(:, ii + 1) + 1 / (beta * h^2) * du; % Incremental correction accelerations

                if sqrt(rr' * rr) / length(rr) < 1.0e-6 % Convergence criteria
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

    %% Function file for the Scanlan nonlinear model
    function [g, ks] = polynomial_nonlinear(u, udot, MM, CC, KK, gamma, beta, h, b, matrixsize, nModes)
        g = (CC - b(5) .* u.^4 - b(4) .* abs(u).^3 - b(3) .* u.^2 - b(2) .* abs(u) - b(1)) * udot + (KK) .* u; % Function value
        kc = CC - b(5) .* u.^4 - b(4) .* abs(u).^3 - b(3) .* u.^2 - b(2) .* abs(u) - b(1);
        kspring = KK + udot .* (-4 .* b(5) .* u^3 - 3 * b(4) * u^3 ./ abs(u) - 2 .* b(3) .* u - b(2) .* u ./ abs(u));
        ks = kspring + gamma .* h ./ (beta .* h.^2) .* kc + 1 ./ (beta .* h.^2) .* MM; % Linearization
    end

    function [g, ks] = limitdampingratio(u, udot, MM, CC, KK, gamma, beta, h, blimit, matrixsize, nModes)
        g = -blimit.*udot+KK*u; % Function value
        g = KK*u; % Function value
        kc = CC - blimit;
        kspring = KK;
        ks = kspring + gamma .* h ./ (beta .* h.^2) .* kc + 1 ./ (beta .* h.^2) .* MM; % Linearization
    end
%     function [g, ks] = limitdampingratioupper(u, udot, MM, CC, KK, gamma, beta, h, blimit, matrixsize, nModes)
%         g = blimit.*udot+KK*u; % Function value
%         g = KK*u; % Function value
%         kc = CC + blimit;
%         kspring = KK;
%         ks = kspring + gamma .* h ./ (beta .* h.^2) .* kc + 1 ./ (beta .* h.^2) .* MM; % Linearization
%     end
end
