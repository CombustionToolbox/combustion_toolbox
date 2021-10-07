function val = entropy_mass(mix)
    % Get the mass specific entropy [kJ/kg-K]
    val = mix.S / mix.mi;
end