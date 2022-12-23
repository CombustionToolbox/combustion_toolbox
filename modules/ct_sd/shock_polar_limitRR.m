function [mix1, mix2, mix2_1, mix3] = shock_polar_limitRR(self, mix1, u1)
    % Obtain polar curves for the given pre-shock conditions
    %
    % Args:
    %    self (struct): Data of the mixture, conditions, and databases
    %    mix1 (struct): Properties of the mixture in the pre-shock state
    %    u1 (float): pre-shock velocity [m/s]
    %
    % Returns:
    %    Tuple containing
    %
    %    * mix1 (struct): Properties of the mixture in the pre-shock state
    %    * mix2 (struct): Properties of the mixture in the post-shock state - polar diagrams from mix1 (incident)
    %    * mix2_1 (struct): Properties of the mixture in the post-shock state - weak shock
    %    * mix3 (struct): Properties of the mixture in the post-shock state - polar diagrams from mix2_1 (reflected)

    % Compute first polar curve
    [mix1, mix2] = shock_polar(self, mix1, u1);

    % Initialization
    mix2_1 = [];
    theta0 = mix2.theta_max / 2; % [deg]
    fprime = -2;
    % Miscellaneous
    it = 0; itMax = self.TN.it_limitRR; STOP = 1.; tol_limitRR = self.TN.tol_limitRR;
    % Loop to find regular reflection limit, which is reached when
    % theta_min_polar2 == 0.
    while STOP > tol_limitRR && it < itMax
        % Update iteration
        it = it + 1;
        % Compute oblique shock for the given deflection angle
        [~, mix2_1] = shock_oblique_theta(self, mix1, u1, theta0, mix2_1);
        % Compute polar curves departing from the post-shock state obtained
        [~, mix3] = shock_polar(self, mix2_1, mix2_1.u);
        % Compute absolute error
        f0 = mix3.theta_max - theta0;
        % Compute first derivative (finite diference - Broyden's method)
        if it > 1
            fprime = (f0 - f0_old) / (theta0 - theta_old);
        end
        % Apply correction
        theta = theta0 - f0 / fprime;
        % Compute STOP criteria
        frel = abs(f0 / theta);
        STOP =  max(abs((theta - theta0) / theta0), frel);
        % Update values
        theta_old = theta0;
        f0_old = f0;
        theta0 = theta;
    end

end