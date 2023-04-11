function self = Miscellaneous()
    % Initialize struct with miscellaneous data
    %
    % Attributes:
    %     description (char): Description of the struct
    %     timer_0 (float): Timer to measure the time of the computations (total time)
    %     timer_loop (float): Timer to measure the time of the computations (only calculation time)
    %     config (struct): Struct with configuration data for plots
    %     FLAG_INITIALIZE (bool): Flag indicating self variable is not fully initialized
    %     FLAG_FIRST (bool): Flag indicating first calculation
    %     FLAG_FOI (bool): % Flag indicating that the reactant mixture has been checked
    %     FLAG_ADDED_SPECIES (bool): Flag indicating that there are added reactants species, because were not considered as products -> to recompute stochiometric matrix
    %     FLAG_N_Fuel (bool): Flag indicating that the number of moles of the fuel species are defined
    %     FLAG_N_Oxidizer (bool): Flag indicating that the number of moles of the oxidant species are defined
    %     FLAG_N_Inert (bool): Flag indicating that the number of moles of the inert species are defined
    %     FLAG_WEIGHT (bool): Flag indicating that the number of moles of the oxidizer/inert speces are defined from its weight percentage
    %     FLAG_RESULTS (bool): Flag to show results in the command window
    %     FLAG_CHECK_INPUTS (bool): Flag indicating that the algorithm has checked the input variables
    %     FLAG_GUI (bool): Flag indicating that the user is using the GUI
    %     FLAG_LABELS (bool): Flag (to be completed)
    %     FLAG_PROP (struct): Struct with flags indicanting if there are several values of the respective property (fieldname)
    %     FLAG_LENGTH (struct): Flag indicating parametric study
    %     export_results (struct): Struct with data to export results
    %     index_LS_original (float): Vector with the indeces original list of products
    %     display_species (struct): Struct with data to display species
    %     i (float): Index of the current calculation
    %
    % Returns:
    %     self (struct): Struct with miscellaneous data
    
    % Description
    self.description = 'Miscellaneous'; 
    % Variables
    % * Timer
    self.timer_0 = [];
    self.timer_loop = [];
    % * Plot
    self.config.position = get_monitor_positions(2); % Default figure position [pixels]
    self.config.innerposition = [0.05 0.05 0.9 0.9]; % Set figure inner position [normalized]
    self.config.outerposition = [0.05 0.05 0.9 0.9]; % Set figure outer position [normalized]
    self.config.linestyle = '-';                     % Set line style for plots
    self.config.symbolstyle = 'o';                   % Set symbol style for plots
    self.config.linewidth = 1.8;                     % Set line width for plots
    self.config.fontsize = 20;                       % Set fontsize
    self.config.colorpalette = 'Seaborn';            % Set Color palette (see brewermap function for more options)
    self.config.colorpaletteLenght = 11;             % Set Maximum number of colors to use in the color palette
    self.config.box = 'off';                         % Display axes outline
    self.config.grid = 'off';                        % Display or hide axes grid lines
    self.config.hold = 'on';                         % Retain current plot when adding new plots
    self.config.axis_x = 'tight';                    % Set x-axis limits
    self.config.axis_y = 'auto';                     % Set y-axis limits
    self.config.xscale = 'linear';                   % Set x-axis scale (linear or log)
    self.config.yscale = 'linear';                   % Set y-axis scale (linear or log)
    self.config.xdir = 'normal';                     % Set x-axis direction (normal or reverse)
    self.config.ydir = 'normal';                     % Set y-axis direction (normal or reverse)
    self.config.title = [];                          % Set title
    self.config.label_type = 'medium';               % Set label with variable (short), name (medium), or name and variable (long)
    self.config.labelx = [];                         % Set x label
    self.config.labely = [];                         % Set y label
    self.config.legend_name = [];                    % Set legend labels
    self.config.legend_location = 'northeastoutside';% Set legend location
    self.config.colorline = [44, 137, 160]/255;      % Default colorline
    self.config.colorlines = [135, 205, 222;...      % Default colorlines
                               95, 188, 211;...      % Default colorlines
                               44, 137, 160;...      % Default colorlines
                               22,  68,  80]/255;    % Default colorlines
    self.config.blue = [0.3725, 0.7373, 0.8275];     % Default colorline
    self.config.gray = [0.50, 0.50, 0.50];           % Default colorline
    self.config.red = [0.64,0.08,0.18];              % Default colorline
    self.config.orange = [212, 85, 0]/255;           % Default colorline
    self.config.brown = [200, 190, 183]/255;         % Default colorline
    self.config.brown2 = [72, 55, 55]/255;           % Default colorline
    % * Flags
    self.FLAG_INITIALIZE = false;                    % Flag indicating self variable is not fully initialized
    self.FLAG_FIRST = true;                          % Flag indicating first calculation
    self.FLAG_FOI = true;                            % Flag indicating that the reactant mixture has to be checked
    self.FLAG_ADDED_SPECIES = false;                 % Flag indicating that there are added reactants species, because were not considered as products -> to recompute stochiometric matrix
    self.FLAG_N_Fuel = true;                         % Flag indicating that the number of moles of the fuel species are defined
    self.FLAG_N_Oxidizer = true;                     % Flag indicating that the number of moles of the oxidant species are defined
    self.FLAG_N_Inert = true;                        % Flag indicating that the number of moles of the inert species are defined
    self.FLAG_WEIGHT = false;                        % Flag indicating that the number of moles of the oxidizer/inert speces are defined from its weight percentage
    self.FLAG_RESULTS = true;                        % Flag to show results in the command window
    self.FLAG_CHECK_INPUTS = false;                  % Flag indicating that the algorithm has checked the input variables
    self.FLAG_GUI = false;                           % Flag indicating that the user is using the GUI
    self.FLAG_LABELS = false;                        % Flag ...
    self.FLAGS_PROP.phi = false;                     % Flag several values of equivalence ratio (phi)
    self.FLAGS_PROP.TR = false;                      % Flag several values of Temperature Reactant (TR)
    self.FLAGS_PROP.TP = false;                      % Flag several values of Temperature Products (TP)
    self.FLAGS_PROP.pR = false;                      % Flag several values of Pressure Reactant (pR)
    self.FLAGS_PROP.pP = false;                      % Flag several values of Pressure Products (pP)
    self.FLAGS_PROP.Aratio = false;                  % Flag several values of Area ratio (Aratio)
    self.FLAGS_PROP.Aratio_c = false;                % Flag several values of Area ratio combustor (Aratio_c)
    self.FLAG_LENGTH_LOOP = false;                   % Flag indicating parametric study
    % * Export data
    self.export_results.value = false;               % Bool variable to export results
    self.export_results.format = '.xls';             % Default fileformat to export results
    self.export_results.filename = 'results';        % Default filename to export results
    % * Others
    self.index_LS_original = [];                     % Indeces original list of products
    self.display_species = {};                       % Species to plot on the molar/mass fraction plots
    self.i = [];                                     % Index of the current computations
end