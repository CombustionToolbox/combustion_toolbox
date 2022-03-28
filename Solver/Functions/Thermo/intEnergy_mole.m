function value = intEnergy_mole(mix)
    % Get the mole specific internal energy [kJ/mol] of the mixture
    %
    % Args:
    %     mix (struct):  Properties of the mixture
    %
    % Returns:
    %     value (float): Mole-basis specific internal energy [kJ/mol] of the mixture

    value = mix.e / mix.N * 1e3;
end