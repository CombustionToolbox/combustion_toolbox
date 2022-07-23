function self = compute_ratio_oxidizers_O2(self)
    % Compute ratio oxidizers/O2
    %
    % Args:
    %     self (struct): Data of the mixture, conditions, and databases
    %
    % Returns:
    %     self (struct): Data of the mixture, conditions, and databases
    
    % Find index O2
    ind_O2 = find_ind(self.PD.S_Oxidizer, 'O2');
    if isempty(ind_O2)
        ind_O2 = find_ind(self.PD.S_Oxidizer, 'O2bLb');
    end
    % Compute proportion inerts/O2
    self.PD.ratio_oxidizers_O2 = self.PD.N_Oxidizer / self.PD.N_Oxidizer(ind_O2);
end