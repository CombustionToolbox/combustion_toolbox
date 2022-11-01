function self = set_react_index(self, species)
    % Set index of no frozen (react) and frozen species
    %
    % Args:
    %     self (struct): Data of the mixture, conditions, and databases
    %     species (str): Frozen species
    %
    % Returns:
    %     self (struct): Data of the mixture, conditions, and databases

    if ~isempty(species)
        index = find_ind(self.S.LS, species);
        self.S.ind_react = [1:index - 1, index + 1:self.S.NS];
        self.S.ind_frozen = index;
    else
        self.S.ind_react = 1:self.S.NS;
        self.S.ind_frozen = [];
    end

end
