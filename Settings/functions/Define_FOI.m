function self = Define_FOI(self, i)
% DEFINE FUEL
if ~self.Misc.FLAG_N_Fuel
    self.PD.N_Fuel = 1;
end
self = Define_F(self);
% DEFINE OXIDIZER
if ~self.Misc.FLAG_N_Oxidizer && ~isempty(self.PD.S_Oxidizer)
    self.PD.N_Oxidizer = self.PD.phi_t/self.PD.phi.value(i);
end
self = Define_O(self);
% DEFINE DILUENTS/INERTS
if ~self.Misc.FLAG_N_Inert && ~isempty(self.PD.S_Inert)
    self.PD.N_Inert = self.PD.phi_t/self.PD.phi.value(i) .* self.PD.proportion_inerts_O2;
end
self = Define_I(self);
% COMPUTE PROPERTIES OF THE REACTIVES FOR THE GIVEN CONDITIONS
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