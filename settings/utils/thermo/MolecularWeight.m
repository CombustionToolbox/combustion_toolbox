function value = MolecularWeight(mix)
    % Get the molecular weight [g/mol] of the mixture
    %
    % Args:
    %     mix (struct):  Properties of the mixture
    %
    % Returns:
    %     value (float): Molecular weight [g/mol] of the mixture

    value = mix.W;
end
