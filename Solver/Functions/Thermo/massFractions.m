function value = massFractions(mix)
    % Get the mass fractions of all the species in the mixture [-]
    %
    % Args:
    %     mix (struct):  Properties of the mixture
    %
    % Returns:
    %     value (float): Mass fractions of all the species in the mixture [-]

    value = mix.Yi;
end