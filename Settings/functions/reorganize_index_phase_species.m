function self = reorganize_index_phase_species(self, LS)
    % Reorginize index of gaseous, condensed and cryogenic species
    self.S.ind_nswt = []; self.S.ind_swt = []; self.S.ind_cryogenic = [];
    self = get_index_phase_species(self, LS);
end