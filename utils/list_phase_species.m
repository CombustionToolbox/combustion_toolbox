function self = list_phase_species(self, LS)
    % Establish cataloged list of species according to the state of the
    % phase (gaseous or condensed). It also obtains the indices of
    % cryogenic liquid species, i.e., liquified gases.
    %
    % Args:
    %     self (struct): Data of the mixture, conditions, and databases
    %     LS (cell):     List of species
    %
    % Returns:
    %     self (struct): Data of the mixture, conditions, and databases

    self = get_index_phase_species(self, LS);
    self.S.LS = self.S.LS([self.S.ind_nswt, self.S.ind_swt]);
    self.S.NS = length(self.S.LS);
    self.S.NG = length(self.S.ind_nswt);
    % Reorginize index of gaseous, condensed and cryogenic species
    self = reorganize_index_phase_species(self, self.S.LS);
end
