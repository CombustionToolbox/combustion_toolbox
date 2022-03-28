function value = enthalpy_mole(mix)
    % Get the mole specific enthalpy [kJ/mol] of the mixture
    %
    % Args:
    %     mix (struct):  Properties of the mixture
    %
    % Returns:
    %     value (float): Mole-basis specific enthalpy [kJ/mol] of the mixture

    value = mix.h / mix.N * 1e3;
end