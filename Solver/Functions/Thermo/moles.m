function val = moles(mix)
    % Get the moles of all the species in the mixture [mol]
    val = mix.Xi * mix.N;
end