function val = intEnergy_mass(mix)
    % Get the mass specific internal energy [kJ/kg]
    val = mix.e / mix.mi;
end