function value = enthalpy_mass(mix)
    % Get the mass specific enthalpy [kJ/kg] of the mixture
    %
    % Args:
    %     mix (struct):  Properties of the mixture
    %
    % Returns:
    %     value (float): Mass-basis specific enthalpy [kJ/kg] of the mixture

    value = mix.h / mix.mi;
end