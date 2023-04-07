function self = set_react_index(self, species)
    % Set index of react (non-frozen) and frozen species
    %
    % Args:
    %     self (struct): Data of the mixture, conditions, and databases
    %     species (char): Frozen species
    %
    % Returns:
    %     self (struct): Data of the mixture, conditions, and databases
    
    % Initialization
    self.S.ind_react = 1:self.S.NS;

    % All species react
    if isempty(species)
        self.S.ind_frozen = [];
        return
    end
    
    % Get index frozen species
    index = find_ind(self.S.LS, species);
    % Set index frozen species
    self.S.ind_frozen = index;
    % Get length initial species
    N_reactants = length([self.PD.S_Fuel, self.PD.S_Oxidizer, self.PD.S_Inert]);
    % Get length frozen species
    N_frozen = length(index);
    % Check if all the species of the reactants are frozen
    if N_frozen == N_reactants
        self.PD.FLAG_FROZEN = true;
        return
    end
    
    for i = 1:N_frozen
        self.S.ind_react(self.S.ind_react == index(i)) = [];
    end

end
