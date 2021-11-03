function self = gui_get_typeSpecies(self)
    typeSpecies = self.UITable_R.Data(:, 4);
    self.ind_Fuel     = contains(typeSpecies, 'Fuel');
    self.ind_Oxidizer = contains(typeSpecies, 'Oxidant');
    self.ind_Inert    = contains(typeSpecies, 'Inert');
end