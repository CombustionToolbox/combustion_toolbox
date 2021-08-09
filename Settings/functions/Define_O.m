function self = Define_O(self)
if ~isempty(self.PD.S_Oxidizer)
    self.PD.R_Oxidizer = SetSpecies(self, self.PD.S_Oxidizer, self.PD.N_Oxidizer, self.PD.TR.value);
else
    self.PD.R_Oxidizer = 0;
end
