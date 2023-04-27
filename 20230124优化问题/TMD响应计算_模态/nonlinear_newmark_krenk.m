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
        u(:, ii + 1) = u(:, ii) + h * udot(:, ii) +1/2 * h ^ 2 * u2dot(:, ii); % Predicted displacements
        Nit = 0; % Number of iterations
        konv = 0; % Has the iterations converged 0=No, 1=Yes
        % Reusidal calculation, system matrices and increment correction

        while Nit < 1000 && konv == 0
            Nit = Nit + 1;

            if Nit > 990
                disp("迭代次数为："+num2str(Nit) + "，已经接近上限")
            end

            [g, Ks] = gfun(u(:, ii + 1), udot(:, ii + 1)); % Calculate function value and the tangent

            % [g, Ks] = gfun(u(:,ii+1),udot(:,ii+1)); % Calculate function value and the tangent
            rr = pp(:, ii + 1) - MM * u2dot(:, ii + 1) - g; % Calculate residual
            du = Ks \ rr; % Increment correction
            u(:, ii + 1) = u(:, ii + 1) + du; % Add incremental correction to the displacement
            udot(:, ii + 1) = udot(:, ii + 1) + gamma * h / (beta * h ^ 2) * du; % Incremental correction for velocities
            u2dot(:, ii + 1) = u2dot(:, ii + 1) + 1 / (beta * h ^ 2) * du; % Incremental correction accelerations

            if sqrt(rr' * rr) / length(rr) < 1.0e-4 % Convergence criteria
                konv = 1; % konv = 1 if the convergence criteria is fulfilled.
            end

        end

    end

end
