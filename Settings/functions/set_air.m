function app = set_air(app, FLAG_IDEAL_AIR)
    % Incluide air in the mixture
    if FLAG_IDEAL_AIR
        app.PD.S_Oxidizer = {'O2'};
        app.PD.S_Inert = {'N2'};
        app.PD.proportion_inerts_O2 = 79/21;
    else
        app.PD.S_Oxidizer = {'O2'};
        app.PD.S_Inert = {'N2', 'Ar', 'CO2'};
        app.PD.proportion_inerts_O2 = [78.084, 0.9365, 0.0319] ./ 20.9476;
    end
end