function value = cv_mass(mix)
    % Get the mass-basis specific heat at constant volume [kJ/kg-K] of the mixture
    %
    % Args:
    %     mix (struct):  Properties of the mixture
    %
    % Returns:
    %     value (float): Mass-basis specific heat at constant volume [kJ/kg-K] of the mixture
    
    value = mix.cV / mix.mi * 1e-3;
end