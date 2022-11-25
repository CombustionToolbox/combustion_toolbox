function self = Species()
    % Initialize struct with problem solution data
    % 
    % Returns:
    %     self (struct): struct with problem solution data

    % Description
    self.description = "Data of the chemical species";
    % Variables
    % * Species and lengths
    self.LS_DB = [];      % List species Database
    self.NS_DB = [];      % Number species Database
    self.NG = [];         % Number gaseous species in the mixture
    self.NS = [];         % Number species in the mixture
    self.LS = [];         % List of species in the mixture
    self.LS_formula = []; % Formula of each species contained in LS
    self.LS_fixed = {'CO2','CO','H2O','H2','O2','N2','He','Ar'}; 
    self.NS_fixed = length(self.LS_fixed); % Number fixed species
    % * Index values
    self.ind_nswt = [];      % Index gaseous species
    self.ind_swt = [];       % Index condensed species
    self.ind_cryogenic = []; % Index cryogenic liquified species
    self.ind_ox_ref = [];    % Index reference oxidizer (default: O2)
    % self.ind_O2 = [];        % Index O2 in LS (deprecated)
    self.ind_fixed = [];     % Index fixed species in LS
    self.ind_all = [];       % Index all species in LS
    self.ind_ions = [];      % Index all ion species in LS
    % * Default list of species for complete combustion
    self.LS_lean = {'CO2', 'H2O', 'N2', 'Ar', 'O2'};       % List of species for a lean complete combustion (equivalence ratio < 1)
    self.LS_rich = {'CO2', 'H2O', 'N2', 'Ar', 'CO', 'H2'}; % List of species for a lean complete combustion (equivalence ratio > 1)
    self.LS_soot = {'N2', 'Ar', 'CO', 'H2', 'Cbgrb', 'CO2', 'H2O'}; % List of species for a lean complete combustion (equivalence ratio > equivalence ratio soot)
    % * Flags
    self.FLAG_COMPLETE = false;
end