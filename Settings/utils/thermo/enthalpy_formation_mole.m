function value = enthalpy_formation_mole(mix)
    % Get the mole specific enthalpy formation [kJ/mol] of the mixture
    %
    % Args:
    %     mix (struct):  Properties of the mixture
    %
    % Returns:
    %     value (float): Mole-basis specific enthalpy formation [kJ/mol] of the mixture

    value = mix.hf / mix.N * 1e3;
end
