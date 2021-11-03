function self = gui_get_reactants(self)
    % Function that obtains the number of moles and the name of each
    % species in the mixture and catalog them in Fuel, Oxidizer or Inert.
    
    % Get indexes of Fuel, Oxidizer and Inert species
    self = gui_get_typeSpecies(self);
    % Fuel
    [self.PD.S_Fuel, self.PD.N_Fuel] = gui_get_reactant(self, self.ind_Fuel);
    % Oxidizer
    [self.PD.S_Oxidizer, self.PD.N_Oxidizer] = gui_get_reactant(self, self.ind_Oxidizer);
    % Inert
    [self.PD.S_Inert, self.PD.N_Inert] = gui_get_reactant(self, self.ind_Inert);
end

% SUB-PASS FUNCTIONS
function [S, N] = gui_get_reactant(self, ind)
    if sum(self.(ind))
        S = regexprep(self.UITable_R.Data(self.(ind), 1),'\s+','')';
        N = self.UITable_R.Data{self.(ind), 2};
    else
        S = [];
        N = [];
    end
end