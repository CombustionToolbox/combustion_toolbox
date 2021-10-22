function self = get_FLAG_N(self)
    % Flag if the number of moles of fuel, oxidant and inert species
    % is specified. If not, consider 1 mole for the fuel and calculate
    % the remaining moles from the equivalence relation.
    if isempty(self.PD.N_Fuel)
        self.Misc.FLAG_N_Fuel = false;
    end
    if isempty(self.PD.N_Oxidizer)
        self.Misc.FLAG_N_Oxidizer = false;
    end
    if isempty(self.PD.N_Inert)
        self.Misc.FLAG_N_Inert = false;
    end
end