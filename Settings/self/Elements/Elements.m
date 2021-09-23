function self = Elements()
   self.description = "Data of the chemical elements";
   self.elements = set_elements();
   self.NE = [];
   self.ind_C = [];
   self.ind_H = [];
   self.ind_O = [];
   self.ind_N = [];
   self.ind_E = []; % Electrons
end