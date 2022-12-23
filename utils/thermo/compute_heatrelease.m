function [q, Q] = compute_heatrelease(mix1, mix2)
    % Compute heat release [J/kg] of the chemical transformation of the mixture 1 to mixture 2
    %
    % Args:
    %     mix1 (struct): Properties of the initial mixture
    %     mix2 (struct): Properties of the final mixture
    %
    % Returns:
    %     q (float): heat release [J/kg] == [m^2/s^2]
    %     Q (float): dimensionless heat release

    
    sound_1 = soundspeed(mix1);
    gamma_1 = adiabaticIndex_sound(mix1);
    mach_1 = velocity_relative(mix1) / sound_1;
    gamma_2 = adiabaticIndex_sound(mix2);

    q = sound_1^2 * ( (1 + gamma_1 * mach_1^2)^2 / (2 * (gamma_2^2 - 1) ) ...
        * (gamma_2 / gamma_1)^2 / mach_1^2 - 1/(gamma_1 - 1) ...
        - mach_1^2 / 2 );
    Q = (gamma_2^2 - 1) / (2 * sound_1^2) * q;
end