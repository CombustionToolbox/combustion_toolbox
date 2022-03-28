function value = intEnergy_mass(mix)
    % Get the mass specific internal energy [kJ/kg] of the mixture
    %
    % Args:
    %     mix (struct):  Properties of the mixture
    %
    % Returns:
    %     value (float): Mass-basis specific internal energy [kJ/kg] of the mixture

    value = mix.e / mix.mi;
end