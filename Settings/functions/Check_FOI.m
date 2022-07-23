function self = Check_FOI(self, FOI_species)
    % Check that fuel species are contained in the list of products (only for initial computations)
    %
    % Args:
    %     self (struct): Data of the mixture, conditions, and databases
    %     FOI_species (bool): Species in the initial mixture (Fuel, Oxidizer, Inert)
    %
    % Returns:
    %     self (struct): Data of the mixture, conditions, and databases

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
        
        % Set index of no frozen (react) and frozen species
        self = set_react_index(self, self.PD.S_Inert);
        
        % Check if O2 is as oxidizer, if not try for O2 in liquid state
        if isempty(find_ind(self.PD.S_Oxidizer, 'O2'))
            self.S.ind_O2 = find_ind(self.S.LS, 'O2bLb');
        end

        self.Misc.FLAG_FIRST = false;
    end
end