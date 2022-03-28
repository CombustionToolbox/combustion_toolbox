function value = entropy_mole(mix)
    % Get the mole specific entropy [kJ/mol-K] of the mixture
    %
    % Args:
    %     mix (struct):  Properties of the mixture
    %
    % Returns:
    %     value (float): Mole-basis specific entropy [kJ/mol-K] of the mixture
    
    value = mix.S / mix.N * 1e3;
end