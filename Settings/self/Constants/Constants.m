function self = Constants()
    self.description = "Constants and tolerances";
    self.R0 = 8.31446261815324;  % [J/(K mol)]. Universal gas constant
    self.A0.description = "Stoichiometric Matrix: number of atoms of each element contained in each species";
    self.A0.value = [];
    self.M0.description = "Matrix with properties of each species";
    self.M0.value = [];
    self.N_prop.description = 'Number of properties in the Matrix properties';
    self.N_prop.value = 12;
    self.N0.description = "Reduced Matrix with number of moles and swtCondensated of each species";
    self.N0.value = [];
    self.MassorMolar = 'mass';
    self.firstrow = true;
    self.mintol_display = 1e-14;
    self.l_phi = []; % length phi vector
end
