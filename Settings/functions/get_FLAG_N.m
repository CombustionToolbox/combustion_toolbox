function self = get_FLAG_N(self)
    % Flag if the number of moles of fuel, oxidant and inert species
    % is specified. If not, consider 1 mole for the fuel and calculate
    % the remaining moles from the equivalence relation.
    self = get_FLAG(self, 'FLAG_N_Fuel', 'N_Fuel');
    self = get_FLAG(self, 'FLAG_N_Oxidizer', 'N_Oxidizer');
    self = get_FLAG(self, 'FLAG_N_Inert', 'N_Inert');
end

% SUB-PASS FUNCTIONS
function self = get_FLAG(self, FLAG_name, cond_name)
    if isempty(self.PD.(cond_name))
        self.Misc.(FLAG_name) = false;
    else
        self.Misc.(FLAG_name) = true;
    end
end