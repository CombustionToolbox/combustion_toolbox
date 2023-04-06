function self = Elements()
    % Initialize struct with elements data
    % 
    % Attributes:
    %     description (char): Description of the struct
    %     elements (cell): Cell with the elements in the periodic table
    %     NE (float): Number of elements
    %     ind_C (float): Index element Carbon
    %     ind_H (float): Index element Hydrogen
    %     ind_O (float): Index element Oxygen
    %     ind_N (float): Index element Nytrogen
    %     ind_E (float): Index element Electron
    %     ind_S (float): Index element Sulfur
    %     ind_Si (float): Index element Silicon
    %
    % Returns:
    %     self (struct): Struct with elements data
    
    % Description
    self.description = 'Data of the chemical elements';
    % Load cell with the elements in the periodic table
    self.elements = set_elements();
    % Variables
    self.NE = [];     % Number of elements
    self.ind_C = [];  % Index element Carbon
    self.ind_H = [];  % Index element Hydrogen
    self.ind_O = [];  % Index element Oxygen
    self.ind_N = [];  % Index element Nytrogen
    self.ind_E = [];  % Index element Electron
    self.ind_S = [];  % Index element Sulfur
    self.ind_Si = []; % Index element Silicon
end