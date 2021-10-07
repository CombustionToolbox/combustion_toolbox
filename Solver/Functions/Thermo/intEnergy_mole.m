function val = intEnergy_mole(mix)
    % Get the mole specific internal energy [kJ/kmol]
    val = mix.e / mix.Ni * 1e3;
end