function self = Define_FOI(self, i)
% COMPUTE PROPERTIES OF THE REACTIVES FOR THE GIVEN CONDITIONS
R = self.PD.R_Fuel + self.PD.R_Oxidizer + self.PD.R_Inert;
self.PS.strR{i} = ComputeProperties(self, R, self.PD.pR.value, self.PD.TR.value);
self.PS.strR{i}.phi = self.PD.phi.value(i);