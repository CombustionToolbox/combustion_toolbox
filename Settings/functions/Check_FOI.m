function self = Check_FOI(self, FOI_species)
    % Check that fuel species are contained in the list of products (only for initial computations)
    if self.Misc.FLAG_FOI
        self.Misc.FLAG_FOI = false;
        self.Misc.LS_original = self.S.LS;
        for i=1:length(FOI_species)
            if ~strcmp(self.S.LS, FOI_species(i))                 
                self.S.LS = [self.S.LS, FOI_species(i)];
                self.S.NS = length(self.S.LS);
                self.S.LS_formula = [self.S.LS_formula, self.DB.(FOI_species{i}).txFormula];
                self.Misc.FLAG_ADDED_SPECIES = true;
            end
        end
        if self.Misc.FLAG_ADDED_SPECIES
            self.S.ind_nswt = []; self.S.ind_swt = []; self.S.ind_cryogenic = [];
            self = ContainedElements(self);
            self = Initialize(self);
            self.Misc.index_LS_original = find_ind(self.S.LS, self.Misc.LS_original);
            % Set indexes phase species to original List Species (for computations)
            self = reorganize_index_phase_species(self, self.Misc.LS_original);
        else
            self.Misc.index_LS_original = 1:1:self.S.NS;
        end
        self.Misc.FLAG_FIRST = false;
    end
end