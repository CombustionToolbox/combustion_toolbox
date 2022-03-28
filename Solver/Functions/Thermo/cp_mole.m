function value = cp_mole(mix)
    % Get the mole-basis specific heat at constant pressure [kJ/mol-K] of the mixture
    %
    % Args:
    %     mix (struct):  Properties of the mixture
    %
    % Returns:
    %     value (float): Mole-basis specific heat at constant pressure [kJ/mol-K] of the mixture
    
    value = mix.cP / mix.N;
end