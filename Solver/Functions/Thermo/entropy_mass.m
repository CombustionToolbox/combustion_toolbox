function value = entropy_mass(mix)
    % Get the mass specific entropy [kJ/kg-K] of the mixture
    %
    % Args:
    %     mix (struct):  Properties of the mixture
    %
    % Returns:
    %     value (float): Mass-basis specific entropy [kJ/kg-K] of the mixture
    
    value = mix.S / mix.mi;
end