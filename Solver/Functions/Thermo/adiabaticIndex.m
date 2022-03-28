function value = adiabaticIndex(mix)
    % Get the adiabatic index [-] of the mixture from the ratio of the specific heat capacities
    %
    % Args:
    %     mix (struct):  Properties of the mixture
    %
    % Returns:
    %     value (float): Adiabatic index of the mixture [-]
    
    value = mix.gamma;
    if isnan(value)
        value = Inf;
    end
end