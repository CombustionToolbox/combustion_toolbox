function self = list_phase_species(self, LS)
    % Establish cataloged list of species according to the state of the 
    % phase (gaseous or condensed). It also obtains the indices of 
    % cryogenic liquid species, e.g., liquified gases.
    self = get_index_phase_species(self, LS);
    if ~isempty(self.S.ind_nswt)
        self.S.ind_nswt = unique(self.S.ind_nswt);
    end
    if ~isempty(self.S.ind_swt)
        self.S.ind_swt  = unique(self.S.ind_swt);
    end
    if ~isempty(self.S.ind_cryogenic)
        self.S.ind_cryogenic = unique(self.S.ind_cryogenic);
    end
    self.S.LS = self.S.LS([self.S.ind_nswt, self.S.ind_swt]);
    self.S.NS = length(self.S.LS);
    self.S.NG = length(self.S.ind_nswt);
    % Reorginize index of gaseous, condensed and cryogenic species
    self = reorganize_index_phase_species(self, self.S.LS);
end