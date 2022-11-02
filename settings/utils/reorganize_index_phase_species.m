function self = reorganize_index_phase_species(self, LS)
    % Reorginize index of gaseous, condensed and cryogenic species
    %
    % Args:
    %     self (struct): Data of the mixture, conditions, and databases
    %     LS (cell): Name list species / list of species
    %
    % Returns:
    %     self (struct):   Data of the mixture, conditions, and databases

    self.S.ind_nswt = []; self.S.ind_swt = []; self.S.ind_cryogenic = [];
    self = get_index_phase_species(self, LS);
end
