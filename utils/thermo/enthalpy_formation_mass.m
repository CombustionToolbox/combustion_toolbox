function value = enthalpy_formation_mass(mix)
    % Get the mass specific enthalpy formation [kJ/kg] of the mixture
    %
    % Args:
    %     mix (struct):  Properties of the mixture
    %
    % Returns:
    %     value (float): Mass-basis specific enthalpy formation [kJ/kg] of the mixture

    value = mix.hf / mix.mi;
end
