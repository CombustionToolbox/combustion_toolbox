function value = moles(mix)
    % Get the moles [mol] of all the species in the mixture
    %
    % Args:
    %     mix (struct):  Properties of the mixture
    %
    % Returns:
    %     value (float): Moles [mol] of all the species in the mixture

    value = mix.Xi * mix.N;
end