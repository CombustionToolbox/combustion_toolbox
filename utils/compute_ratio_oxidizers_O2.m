function self = compute_ratio_oxidizers_O2(self)
    % Compute ratio oxidizers/O2
    %
    % Args:
    %     self (struct): Data of the mixture, conditions, and databases
    %
    % Returns:
    %     self (struct): Data of the mixture, conditions, and databases

    % Compute proportion inerts/O2
    self.PD.ratio_oxidizers_O2 = self.PD.N_Oxidizer / self.PD.N_Oxidizer(self.S.ind_ox_ref);
end
