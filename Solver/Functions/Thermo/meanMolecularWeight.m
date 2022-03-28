function value = meanMolecularWeight(mix)
    % Get the mean molecular weight [g/mol] of the mixture
    %
    % Args:
    %     mix (struct):  Properties of the mixture
    %
    % Returns:
    %     value (float): Mean molecular weight [g/mol] of the mixture

    value = mix.W;
end