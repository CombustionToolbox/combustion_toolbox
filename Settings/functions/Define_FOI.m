function self = Define_FOI(self, i)
    % Set up mixture: fuel, oxidizer and diluent/inert species

    % Check reactant species are contained in the list of products (initial computations)
    self = Check_FOI(self, [self.PD.S_Fuel, self.PD.S_Oxidizer, self.PD.S_Inert]);
    % Define Fuel
    if ~self.Misc.FLAG_N_Fuel
        self.PD.N_Fuel = 1;
    end
    self = Define_F(self);
    % Define Oxidizer
    if ~self.Misc.FLAG_N_Oxidizer && ~isempty(self.PD.S_Oxidizer)
        self.PD.N_Oxidizer = self.PD.phi_t/self.PD.phi.value(i);
    end
    self = Define_O(self);
    % Define Diluent/Inert
    if ~self.Misc.FLAG_N_Inert && ~isempty(self.PD.S_Inert)
        self.PD.N_Inert = self.PD.phi_t/self.PD.phi.value(i) .* self.PD.proportion_inerts_O2;
    end
    self = Define_I(self);
    % Compute properties of the reactives for the given conditions
    R = self.PD.R_Fuel + self.PD.R_Oxidizer + self.PD.R_Inert;
    self.PS.strR{i} = ComputeProperties(self, R, self.PD.pR.value, self.PD.TR.value);
    self.PS.strR{i}.phi = self.PD.phi.value(i);
    self.PS.strR{i}.LS  = merged_cells({self.PD.S_Fuel, self.PD.S_Oxidizer, self.PD.S_Inert});
    [~, ind_LS, ~] = intersect(self.PS.strR{i}.LS, self.S.LS);
    self.PS.strR{i}.ind_LS = ind_LS;
end

% SUB-PASS FUNCTIONS
function merged = merged_cells(cells)
    merged = [];
    for i = 1:length(cells)
        merged = [merged, cells{i}];
    end
end