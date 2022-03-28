function value = gibbs_mole(mix)
    % Get the mole specific gibbs free energy [kJ/mol] of the mixture
    %
    % Args:
    %     mix (struct):  Properties of the mixture
    %
    % Returns:
    %     value (float): Mole-basis specific gibbs free energy [kJ/mol] of the mixture

    value = mix.g / mix.N * 1e3;
end