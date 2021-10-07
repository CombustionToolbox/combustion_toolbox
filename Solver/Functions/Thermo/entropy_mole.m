function val = entropy_mole(mix)
    % Get the mole specific entropy [kJ/kmol-K]
    val = mix.S / mix.Ni * 1e3;
end