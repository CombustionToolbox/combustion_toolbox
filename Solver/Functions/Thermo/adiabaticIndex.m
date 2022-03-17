function val = adiabaticIndex(mix)
    % Get the adiabatic index [-]
    val = mix.gamma;
    if isnan(val)
        val = Inf;
    end
end