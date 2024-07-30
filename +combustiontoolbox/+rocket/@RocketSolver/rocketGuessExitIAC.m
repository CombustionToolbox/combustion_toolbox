function log_Pe = rocketGuessExitIAC(mix2, mix3, areaRatio, FLAG_SUBSONIC)
    % Compute guess logarithm of the ratio pressure_inf / pressure_exit
    % [-] for the given Area ratio [-] and indicanting if the point of
    % interest is in the subsonic area ratios or the supersonic area ratios
    %
    % This method is based on the method outlined in Gordon, S., & McBride,
    % B. J. (1994). NASA reference publication, 1311.
    %
    % Args:
    %      mix2 (Mixture): Properties of the mixture at the outlet of the chamber
    %      mix3 (Mixture): Properties of the mixture at the throat
    %      areaRatio (float): Ratio area_exit / area_throat
    %      FLAG_SUBSONIC (bool): Flag indicating if the areaRatio refer to the subsonic region or the supersonic region
    %
    % Returns:
    %      log_P (float): Log pressure ratio [-]
    %
    % Example:
    %      log_P = rocketGuessExitIAC(mix2, mix3, 3, false)

    % Definitions
    log_Pt = log(mix2.p / mix3.p);

    % Compute guess
    if FLAG_SUBSONIC
        log_Pe = computePressureSubsonic(log_Pt, areaRatio);
    else
        log_Pe = computePressureSupersonic(log_Pt, areaRatio, mix2.gamma_s);
    end

end

% SUB-PASS FUNCTIONS
function log_Pe = computePressureSubsonic(log_Pt, areaRatio)
    % Compute guess logarithm of the pressure ratio for the subsonic area
    % ratios
    if areaRatio >= 1.09
        log_Pe = log_Pt / (areaRatio + 10.587 * log(areaRatio)^3 + 9.454 * log(areaRatio));
    elseif areaRatio > 1.0001
        log_Pe = 0.9 * log_Pt / (areaRatio + 10.587 * log(areaRatio)^3 + 9.454 * log(areaRatio));
    else
        error('The area ratio can not be smaller or equal than 1.0001');
    end

end

function log_Pe = computePressureSupersonic(log_Pt, areaRatio, gamma_s)
    % Compute guess logarithm of the pressure ratio for the subsonic area
    % ratios
    if areaRatio >= 2
        log_Pe = gamma_s + 1.4 * log(areaRatio);
    elseif areaRatio > 1.0001
        log_Pe = log_Pt + sqrt(3.294 * (log(areaRatio)^2) + 1.535 * log(areaRatio));
    else
        error('The area ratio can not be smaller or equal than 1.0001');
    end

end