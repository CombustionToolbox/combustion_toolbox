function self = Species()
    % Initialize struct with chemical species data
    %
    % Attributes:
    %     description (char): Description of the struct
    %     LS_DB (cell): List of species in the database
    %     NS_DB (float): Number of species in the database
    %     NG (float): Number of gaseous species in the mixture
    %     NS (float): Number of species in the mixture
    %     LS (cell): List of species in the mixture
    %     LS_formula (cell): Formula of each species contained in LS
    %     ind_nswt (float): Indeces gaseous species
    %     ind_swt (float): Indeces condensed species
    %     ind_cryogenic (float): Indeces cryogenic liquified species
    %     ind_ox_ref (float): Indeces reference oxidizer (default: O2)
    %     ind_ions (float): Indeces ionized species in LS
    %     ind_react (float): Indeces react species
    %     ind_frozen (float): Indeces inert/frozen species
    %     LS_lean (cell): List of species for a lean complete combustion (equivalence ratio < 1)
    %     LS_rich (cell): List of species for a lean complete combustion (equivalence ratio > 1)
    %     LS_soot (cell): List of species for a lean complete combustion (equivalence ratio > equivalence ratio soot)
    %     FLAG_COMPLETE (bool): Flag indicating if the complete combustion is considered
    %     FLAG_BURCAT (bool): Find all the combinations of species from the database (without BURCAT's DB) that can appear as products for the given list of reactants
    %     FLAG_ION (bool): Flag indicating to include ionized species in the automatic finder of species
    %
    % Returns:
    %     self (struct): struct with chemical species data

    % Description
    self.description = 'Data of the chemical species';
    % Variables
    % * Species and lengths
    self.LS_DB = [];      % List species Database
    self.NS_DB = [];      % Number species Database
    self.NG = [];         % Number gaseous species in the mixture
    self.NS = [];         % Number species in the mixture
    self.LS = [];         % List of species in the mixture
    self.LS_formula = []; % Formula of each species contained in LS
    % * Index values
    self.ind_nswt = [];      % Indeces gaseous species
    self.ind_swt = [];       % Indeces condensed species
    self.ind_cryogenic = []; % Indeces cryogenic liquified species
    self.ind_ox_ref = [];    % Indeces reference oxidizer (default: O2)
    % self.ind_O2 = [];        % Indeces O2 in LS (deprecated)
    % self.ind_fixed = [];     % Indeces fixed species in LS (deprecated)
    % self.ind_all = [];       % Indeces all species in LS (deprecated)
    self.ind_ions = [];      % Indeces ionized species in LS
    self.ind_react = [];     % Indeces react species
    self.ind_frozen = [];    % Indeces inert/frozen species
    % * Default list of species for complete combustion
    self.LS_lean = {'CO2', 'H2O', 'N2', 'Ar', 'O2'};       % List of species for a lean complete combustion (equivalence ratio < 1)
    self.LS_rich = {'CO2', 'H2O', 'N2', 'Ar', 'CO', 'H2'}; % List of species for a lean complete combustion (equivalence ratio > 1)
    self.LS_soot = {'N2', 'Ar', 'CO', 'H2', 'Cbgrb', 'CO2', 'H2O'}; % List of species for a lean complete combustion (equivalence ratio > equivalence ratio soot)
    % * Flags
    self.FLAG_COMPLETE = false;
    self.FLAG_BURCAT = false; % Find all the combinations of species from the database (without BURCAT's DB) that can appear as products for the given list of reactants
    self.FLAG_ION = false; % Flag indicating to include ionized species in the automatic finder of species
end