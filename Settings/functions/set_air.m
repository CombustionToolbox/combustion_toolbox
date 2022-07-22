function self = set_air(self, FLAG_IDEAL_AIR)
    % Incluide air in the mixture
    if FLAG_IDEAL_AIR
        self.PD.S_Oxidizer = {'N2', 'O2'};
        self.PD.ratio_oxidizers_O2 = [79, 21]/21;
    else
        self.PD.S_Oxidizer = {'N2', 'O2', 'Ar', 'CO2'};
        self.PD.ratio_oxidizers_O2 = [78.084, 20.9476, 0.9365, 0.0319] ./ 20.9476;
    end
end