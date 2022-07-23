function self = Define_O(self)
    % Set Oxidizer of the mixture
    %
    % Args:
    %     self (struct): Data of the mixture, conditions, and databases
    %
    % Returns:
    %     self (struct): Data of the mixture, conditions, and databases
    
    if ~isempty(self.PD.S_Oxidizer)
        self.PD.R_Oxidizer = SetSpecies(self, self.PD.S_Oxidizer, self.PD.N_Oxidizer, self.PD.TR.value);
    else
        self.PD.R_Oxidizer = 0;
    end
end