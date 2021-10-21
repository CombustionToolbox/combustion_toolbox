function self = Define_FOI(self, i)
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
    merged = cell(1, numel(cells));
    for j = 1:length(cells)
        for i = 1:length(cells{j})
            merged(j+i-1) = cells{i, j};
        end
    end
end