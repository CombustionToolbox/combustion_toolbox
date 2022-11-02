function self = get_index_phase_species(self, LS)
    % Get index of gaseous, condensed and cryogenic species
    %
    % Args:
    %     self (struct): Data of the mixture, conditions, and databases
    %     LS (cell): Name list species / list of species
    %
    % Returns:
    %     self (struct):   Data of the mixture, conditions, and databases

    % Initialize vectors
    self.S.ind_nswt = [];
    self.S.ind_swt = [];
    self.S.ind_cryogenic = [];
    self.S.ind_ions = [];

    % Get index
    for ind = 1:length(LS)
        Species = LS{ind};

        if ~self.DB.(Species).phase
            self.S.ind_nswt = [self.S.ind_nswt, ind];
        else
            self.S.ind_swt = [self.S.ind_swt, ind];

            if ~self.DB.(Species).ctTInt
                self.S.ind_cryogenic = [self.S.ind_cryogenic, ind];
            end

        end

    end

    self.S.ind_ions = get_index_ions(self.S.LS);
end
