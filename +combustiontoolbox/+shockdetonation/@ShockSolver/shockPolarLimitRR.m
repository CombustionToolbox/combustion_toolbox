function [mix1, mix2, mix2_1, mix3] = shockPolarLimitRR(obj, mix1, u1)
    % Obtain polar curves for the given pre-shock conditions using
    % Broyden's method
    %
    % Args:
    %     obj (ShockSolver): ShockSolver object
    %     mix1 (Mixture):    Properties of the mixture in the pre-shock state
    %     u1 (float):        Pre-shock velocity [m/s]
    %
    % Returns:
    %    Tuple containing
    %
    %    * mix1 (Mixture): Properties of the mixture in the pre-shock state
    %    * mix2 (Mixture): Properties of the mixture in the post-shock state - polar diagrams from mix1 (incident)
    %    * mix2_1 (Mixture): Properties of the mixture in the post-shock state - weak shock
    %    * mix3 (Mixture): Properties of the mixture in the post-shock state - polar diagrams from mix2_1 (reflected)
    %
    % Example:
    %    [mix1, mix2, mix2_1, mix3] = shockPolarLimitRR(ShockSolver(), mix1, u1)
    
    % Miscellaneous
    temp = obj.FLAG_RESULTS;
    obj.FLAG_RESULTS = false;

    % Compute first polar curve
    [mix1, mix2] = shockPolar(obj, mix1, u1);

    % Initialization
    mix2_1 = [];
    theta0 = mix2.thetaMax / 2; % [deg]
    fprime = -2;

    % Miscellaneous
    it = 0; itMax = obj.itLimitRR; STOP = 1; tolLimitRR = obj.tolLimitRR;
    
    % Create a temporal shallow copy of mix1 to avoid overwrite results
    temp_mix1 = mix1.copy;

    % Loop to find regular reflection limit, which is reached when
    % theta_min_polar2 == 0
    while STOP > tolLimitRR && it < itMax
        % Update iteration
        it = it + 1;
        % Compute oblique shock for the given deflection angle
        [~, mix2_1] = shockObliqueTheta(obj, temp_mix1, u1, theta0, mix2_1);
        % Compute polar curves departing from the post-shock state obtained
        [~, mix3] = shockPolar(obj, mix2_1, mix2_1.u);
        % Compute absolute error
        f0 = mix3.thetaMax - theta0;
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

    % Miscellaneous
    obj.FLAG_RESULTS = temp;
end