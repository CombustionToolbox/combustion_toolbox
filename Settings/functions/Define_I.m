function self = Define_I(self)
if ~isempty(self.PD.S_Inert)
    self.PD.R_Inert = SetSpecies(self, self.PD.S_Inert, self.PD.N_Inert, self.PD.TR.value);
else
    self.PD.R_Inert = 0; % case without inert gases
end
if ~isempty(self.PD.S_Oxidizer) || ~isempty(self.PD.S_Inert)
    self.PS.strR_Oxidizer = ComputeProperties(self, self.PD.R_Oxidizer + self.PD.R_Inert, self.PD.pR.value, self.PD.TR.value);
end