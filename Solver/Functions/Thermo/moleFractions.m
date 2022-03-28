function value = moleFractions(mix)
    % Get the mole fractions of all the species in the mixture [-]
    %
    % Args:
    %     mix (struct):  Properties of the mixture
    %
    % Returns:
    %     value (float): Mole fractions of all the species in the mixture [-]

    value = mix.Xi;
end