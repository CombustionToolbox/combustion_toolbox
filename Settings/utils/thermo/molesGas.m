function val = molesGas(mix)
    % Get the moles of the gases in the mixture [mol]
    val = sum(mix.N * mix.Xi(mix.phase == 0));
end
