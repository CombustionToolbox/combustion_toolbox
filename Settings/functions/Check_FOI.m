function self = Check_FOI(self, FOI_species)
    % Check that fuel species are contained in the list of products (only for initial computations)
    if self.Misc.FLAG_FOI
        self.Misc.FLAG_FOI = false;
        self.Misc.index_LS_original = 1:1:self.S.NS;
        for i=1:length(FOI_species)
            if ~strcmp(self.S.LS, FOI_species(i))                 
                self.S.LS = [self.S.LS, FOI_species(i)];
                self.S.NS = length(self.S.LS);
                self.S.LS_formula = [self.S.LS_formula, self.DB.(FOI_species{i}).txFormula];
                self.Misc.FLAG_ADDED_SPECIES = true;
            end
        end
        if self.Misc.FLAG_ADDED_SPECIES
            self.S.ind_nswt = [];
            self.S.ind_swt = [];
            self = ContainedElements(self);
            self = Initialize(self);
            self.S.ind_nswt = intersect(self.S.ind_nswt, self.Misc.index_LS_original);
            self.S.ind_swt = intersect(self.S.ind_swt, self.Misc.index_LS_original);
            self.S.ind_cryogenic = intersect(self.S.ind_cryogenic, self.Misc.index_LS_original);
        end
        self.Misc.FLAG_FIRST = false;
    end
end