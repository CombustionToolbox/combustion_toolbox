function value = cp_mass(mix)
    % Get the mass-basis specific heat at constant pressure [kJ/kg-K] of the mixture
    %
    % Args:
    %     mix (struct):  Properties of the mixture
    %
    % Returns:
    %     value (float): Mass-basis specific heat at constant pressure [kJ/kg-K] of the mixture
    
    value = mix.cP / mix.mi * 1e-3;
end