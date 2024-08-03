function value = entropy_mole(mix)
    % Get the mole specific entropy [kJ/mol-K] of the mixture
    %
    % Args:
    %     mix (struct):  Properties of the mixture
    %
    % Returns:
    %     value (float): Mole-basis specific entropy [kJ/mol-K] of the mixture

    value = mix.s / mix.N * 1e-3;
end
