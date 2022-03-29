function self = Elements()
    % Initialize struct with elements data
    % 
    % Returns:
    %     self (struct): struct with elements data
    
    % Description
    self.description = "Data of the chemical elements";
    % Load cell with the elements in the periodic table
    self.elements = set_elements();
    % Variables
    self.NE = [];    % Number of elements
    self.ind_C = []; % Index element Carbon
    self.ind_H = []; % Index element Hydrogen
    self.ind_O = []; % Index element Oxygen
    self.ind_N = []; % Index element Nytrogen
    self.ind_E = []; % Index element Electron
end