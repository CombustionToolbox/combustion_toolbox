function self = check_FOI(self, FOI_species)
    % Check that fuel species are contained in the list of products (only for initial computations)
    %
    % Args:
    %     self (struct): Data of the mixture, conditions, and databases
    %     FOI_species (bool): Species in the initial mixture (Fuel, Oxidizer, Inert)
    %
    % Returns:
    %     self (struct): Data of the mixture, conditions, and databases

    if ~self.Misc.FLAG_FOI
        return
    end
    
    self.Misc.FLAG_FOI = false;
    self.Misc.LS_original = self.S.LS;

    for i = 1:length(FOI_species)

        if ~strcmp(self.S.LS, FOI_species(i))
            self.S.LS = [self.S.LS, FOI_species(i)];
            self.S.NS = length(self.S.LS);
            self.S.LS_formula = [self.S.LS_formula, self.DB.(FOI_species{i}).txFormula];
            self.Misc.FLAG_ADDED_SPECIES = true;
        end

    end

    if self.Misc.FLAG_ADDED_SPECIES
        self.S.ind_nswt = []; self.S.ind_swt = []; self.S.ind_cryogenic = [];
        self = contained_elements(self);
        self = initialize(self);
        self.Misc.index_LS_original = find_ind(self.S.LS, self.Misc.LS_original);
        % Set indexes phase species to original List Species (for computations)
        self = reorganize_index_phase_species(self, self.Misc.LS_original);
        % Fill constant values of properties matrix
        self = set_species_initilize(self, self.S.LS);
    else
        self.Misc.index_LS_original = 1:1:self.S.NS;
    end

    % Set index of no frozen (react) and frozen species
    self = set_react_index(self, self.PD.S_Inert);

    % Get oxidizer of reference
    self = get_oxidizer_reference(self);

    % Check if oxidizer/inert inputs comes from weight percentage
    if ~isempty(self.PD.wt_ratio_oxidizers) || ~isempty(self.PD.wt_ratio_inerts)
        self.Misc.FLAG_WEIGHT = true;
    end
    
    self.Misc.FLAG_FIRST = false;
end
