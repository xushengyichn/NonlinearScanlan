%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: Shengyi xushengyichn@outlook.com
%Date: 2022-05-15 16:16:48
%LastEditors: Shengyi xushengyichn@outlook.com
%LastEditTime: 2022-06-07 00:42:49
%FilePath: \NonlinearScanlan\NonlinearScanlan_withTMD_Mohamed1996.m
%Description: 本函数目的为根据参考文献，验证nonlinear Scanlan模型施加tmd计算响应的正确性
% Control by passive TMD of wind-induced nonlinear vibrations in cable stayed bridges
% M. Abdel-Rohman and H. Askar
% Journal of Vibration and Control 1996 Vol. 2 Issue 2 Pages 251-267
%
%Copyright (c) 2022 by Shengyi xushengyichn@outlook.com, All Rights Reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a function
close all
clc
clear
% function out = NonlinearScanlan(Fre, Mass, Zeta0, rho, D, U, Y1k, epsilonk, Y2k, ClK, t, P, u0, udot0, mtmd, ctmd, ktmd)

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
% mtmd     = mass of tmd
% ctmd     = damping of tmd
% ktmd     = stiffness of tmd

%--------------------------------------------------------------------------
% beta = 0,     gamma = 1/2 => explicit central difference method
% beta = 1/4,   gamma = 1/2 => undamped trapezoidal rule (implicit)

%--------------------------------------------------------------------------
% Default inspection data
% Marra, A. M., et al. (2015). "Measurements and improved model of
% vortex-induced vibration for an elongated rectangular cylinder."
% Journal of Wind Engineering and Industrial Aerodynamics 147: 358-367.

B = 33/2;
rho = 1.225;
M1 = 1;
H = 14.06 / M1;
epsilon = 237000;
phi2 = 0.12268 * B^2 / epsilon / 1.5376;
phi = sqrt(phi2);

U = 1.5376 / rho / B / H / phi2;
omega = 1.99;

gamma = 1/2; % Factor in the Newmark algorithm
beta = 1/4; % Factor in the Newmark algorithm

casenum = 1;

Fre = omega / 2 / pi;
Mass = M1;

Zeta0 = 0.01;

Y1k = H;

epsilonk = epsilon;

Y2k = 0;

ClK = 0;

h = 0.01;
T = 150;
t = 0:h:T;

P = zeros(2, length(t));
pp = P;

matrixsize = 1; %calculate one degree of freedom

u0 = [0.1; 0];
udot0 = [0; 0];
mtmd = Mass * 0.001;
disp("质量比为：" + num2str(mtmd / Mass * 100) + "%")
Ftmd = Fre*0.98;
Omegatmd = Ftmd * 2 * pi;
ktmd = Omegatmd^2 * mtmd;
ctmd = 0.15 * 2 * mtmd * Omegatmd;
Y1 = Y1k;
epsilon = epsilonk;
Y2 = Y2k;

MM = [Mass 0; 0 mtmd];
CC = zeros(size(MM, 1), size(MM, 2));
CC(1, 1) = Zeta0 * 4 * pi * Mass * Fre + ctmd;
CC(2, 1) = -ctmd;
CC(1, 2) = -ctmd;
CC(2, 2) = ctmd;
KK = zeros(size(MM, 1), size(MM, 2));
KK(1, 1) = 4 * pi^2 * Mass * Fre^2 +ktmd;
KK(2, 1) = -ktmd;
KK(1, 2) = -ktmd;
KK(2, 2) = ktmd;

%% Calculate the response

nModes = 1;

gfun = @(u, udot) Scanlan_nonlinear(u, udot, MM, CC, KK, gamma, beta, h, rho, U, B, Y1, epsilon, Y2,phi);
[u, udot, u2dot] = nonlinear_newmark_krenk(gfun, MM, pp, u0, udot0, gamma, beta, h);

% u_nounit = u / D;
% s = t * U / D;
out = [t; u];
% subplot(1, 2, 1)
% plot(out(1, :), out(2, :))
% % ylim([-0.01 0.01])
% title("main structure")
% subplot(1, 2, 2)
% plot(out(1, :), out(3, :))
% title("TMD")
% % ylim([-0.1 0.1])
% %     figure
% %     plot(s,u_nounit)
% 
% %     figure
% %     [psd_avg, f, psd_plot] = fft_transfer(1/h,u_nounit');
% %     plot(f,psd_plot)

figure
subplot(1, 2, 1)
plot(out(1, :), out(2, :))
xlabel("Time(sec)")
ylabel("Displacement(m)")
out2 = [u; udot];
subplot(1, 2, 2)
plot(out2(3, :), out2(1, :))
ylim([-10 10])
xlim([-20 20])
xlabel("Velocity(m/sec)")
ylabel("Displacement(m)")

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

        end

    end

end

%% Function file for the Scanlan nonlinear model
function [g, ks] = Scanlan_nonlinear(u, udot, MM, CC, KK, gamma, beta, h, rho, U, B, Y1, epsilon, Y2,phi)
    Caero = [-rho*U*B*Y1*phi^2+rho*U/B*Y1*epsilon*phi^4*u(1,1)^2 0; 0 0];

    g = (CC + Caero) * udot + (KK) * u; % Function value
    kc = CC + Caero;

    kspring = KK + [2*rho*U/B*Y1*epsilon*phi^4*u(1,1)*udot(1,1) 0; 0 0];
    ks = kspring + gamma .* h ./ (beta .* h.^2) .* kc + 1 ./ (beta .* h.^2) .* MM; % Linearization
end

% end
