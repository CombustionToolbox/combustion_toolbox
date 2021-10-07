function val = enthalpy_mass(mix)
    % Get the mass specific enthalpy [kJ/kg]
    val = mix.h / mix.mi;
end