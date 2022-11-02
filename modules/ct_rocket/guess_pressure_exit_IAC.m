function log_Pe = guess_pressure_exit_IAC(mix2, mix3, Aratio, FLAG_SUBSONIC)
    % Compute guess logarithm of the ratio pressure_inf / pressure_exit
    % [-] for the given Area ratio [-] and indicanting if the point of
    % interest is in the subsonic area ratios or the supersonic area ratios
    %
    % Args:
    %      mix2 (struct): Properties of the mixture at the outlet of the chamber
    %      mix3 (struct): Properties of the mixture at the throat
    %      Aratio (struct): Ratio area_exit / area_throat
    %      FLAG_SUBSONIC (bool): Flag indicating if the Aratio refer to the subsonic region or the supersonic region
    %
    % Returns:
    %      log_P (float): Log pressure ratio [-]

    % Definitions
    log_Pt = log(mix2.p / mix3.p);

    % Compute guess
    if FLAG_SUBSONIC
        log_Pe = compute_pressure_subsonic(log_Pt, Aratio);
    else
        log_Pe = compute_pressure_supersonic(log_Pt, Aratio, mix2.gamma_s);
    end

end

% SUB-PASS FUNCTIONS
function log_Pe = compute_pressure_subsonic(log_Pt, Aratio)
    % Compute guess logarithm of the pressure ratio for the subsonic area
    % ratios
    if Aratio >= 1.09
        log_Pe = log_Pt / (Aratio + 10.587 * log(Aratio)^3 + 9.454 * log(Aratio));
    elseif Aratio > 1.0001
        log_Pe = 0.9 * log_Pt / (Aratio + 10.587 * log(Aratio)^3 + 9.454 * log(Aratio));
    else
        error('The area ratio can not be smaller or equal than 1.0001');
    end

end

function log_Pe = compute_pressure_supersonic(log_Pt, Aratio, gamma_s)
    % Compute guess logarithm of the pressure ratio for the subsonic area
    % ratios
    if Aratio >= 2
        log_Pe = gamma_s + 1.4 * log(Aratio);
    elseif Aratio > 1.0001
        log_Pe = log_Pt + sqrt(3.294 * (log(Aratio)^2) + 1.535 * log(Aratio));
    else
        error('The area ratio can not be smaller or equal than 1.0001');
    end

end
