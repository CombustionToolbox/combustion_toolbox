function self = compute_proportion_inerts_O2(self)
    % Compute proportion inerts/O2
    
    % Find index O2
    ind_O2 = find_ind(self.PD.S_Oxidizer, 'O2');
    % Compute proportion inerts/O2
    self.PD.proportion_inerts_O2 = self.PD.N_Inert./(sum(self.PD.N_Inert) + self.PD.N_Oxidizer(ind_O2));
end
