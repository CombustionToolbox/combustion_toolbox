function self = define_I(self)
    % Set Inert of the mixture
    %
    % Args:
    %     self (struct): Data of the mixture, conditions, and databases
    %
    % Returns:
    %     self (struct): Data of the mixture, conditions, and databases

    if ~isempty(self.PD.S_Inert)
        self.PD.R_Inert = set_species(self, self.PD.S_Inert, self.PD.N_Inert, self.PD.TR.value);
    else
        self.PD.R_Inert = 0; % case without inert gases
    end

    if ~isempty(self.PD.S_Oxidizer) || ~isempty(self.PD.S_Inert)
        self.PS.strR_Oxidizer = compute_properties(self, self.PD.R_Oxidizer + self.PD.R_Inert, self.PD.pR.value, self.PD.TR.value);
    end

end