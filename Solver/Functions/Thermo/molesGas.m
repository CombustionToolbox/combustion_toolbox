function val = molesGas(mix)
    % Get the moles of the gases in the mixture [mol]
    val = sum(mix.N * mix.Xi(mix.swtCond == 0));
end