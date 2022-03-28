function value = gibbs_mass(mix)
    % Get the mass specific gibbs free energy [kJ/kg] of the mixture
    %
    % Args:
    %     mix (struct):  Properties of the mixture
    %
    % Returns:
    %     value (float): Mass-basis specific gibbs free energy [kJ/kg] of the mixture

    value = mix.g / mix.mi;
end