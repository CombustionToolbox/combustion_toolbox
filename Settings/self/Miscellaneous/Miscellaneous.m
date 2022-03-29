function self = Miscellaneous()
    % Initialize struct with miscellaneous data
    % 
    % Returns:
    %     self (struct): struct with miscellaneous data
    
    % Description
    self.description = "Constants"; 
    % Variables
    %   * Timer
    self.timer_0 = [];
    self.timer_loop = [];
    %   * Plot
    self.config.linewidth = 1.8;              % Default linewidth for plots
    self.config.fontsize = 22;                % Default fontsize
    self.config.colorpalette = 'Spectral';    % Default color palette (see brewermap function for more options)
    self.config.colorpaletteLenght = 12;      % Default number of colors to use in the color palette
    %   * Flags
    self.FLAG_FIRST = true;                   % Flag indicating first calculation
    self.FLAG_FOI = true;                     % Flag indicating that the computations of the reactant mixture ....
    self.FLAG_ADDED_SPECIES = false;          % Flag indicating that there are added reactants species, because were not considered as products -> to recompute stochiometric matrix
    self.FLAG_N_Fuel = true;                  % Flag indicating that the number of moles of the fuel species are defined
    self.FLAG_N_Oxidizer = true;              % Flag indicating that the number of moles of the oxidant species are defined
    self.FLAG_N_Inert = true;                 % Flag indicating that the number of moles of the inert species are defined
    self.FLAG_RESULTS = true;                 % Flag to show results in the command window
    self.FLAG_CHECK_INPUTS = false;           % Flag indicating that the algorithm has checked the input variables
    self.FLAG_GUI = false;                    % Flag indicating that the user is using the GUI
    self.FLAG_LABELS = false;                 % Flag ...
    self.FLAG_PROP = [];                      % Struct with flags indicanting if there are several values of the respective property (fieldname)
    self.FLAGS_PROP.TR = false;               % Flag several values of Temperature Reactant (TR)
    self.FLAG_LENGTH_LOOP = false;            % Flag for several computations
    %   * Export data
    self.export_results.value = false;        % Bool variable to export results
    self.export_results.format = '.xls';      % Default fileformat to export results
    self.export_results.filename = 'results'; % Default filename to export results
    %   * Others
    self.index_LS_original = [];              % Indexes original list of products
    self.display_species = {};                % Species to plot on the molar/mass fraction plots
    self.i = [];                              % Index of the current computations
end