function self = Define_F(self)
    % Set Fuel of the mixture
    if ~isempty(self.PD.S_Fuel)
        self = Check_FOI(self, self.PD.S_Fuel);
        self.PD.R_Fuel = SetSpecies(self, self.PD.S_Fuel, self.PD.N_Fuel, self.PD.TR.value);
        self.PS.strR_Fuel = ComputeProperties(self, self.PD.R_Fuel, self.PD.pR.value, self.PD.TR.value);
        self.PD.Fuel.x = self.PS.strR_Fuel.x;
        self.PD.Fuel.y = self.PS.strR_Fuel.y;
        self.PD.Fuel.z = self.PS.strR_Fuel.z;
        self.PD.Fuel.w = self.PS.strR_Fuel.w;
        
        self.PD.phi_t = self.PD.Fuel.x+self.PD.Fuel.y/4-self.PD.Fuel.z/2;
    else
        self.PD.R_Fuel = 0; self.PD.phi_t = 1;
        self.PD.Fuel.x = 0; self.PD.Fuel.y = 0;
        self.PD.Fuel.z = 0; self.PD.Fuel.w = 0;
        self.C.FLAG_Fuel = 0;
    end
end