function self = App(varargin)
    % Generate self variable with all the data required to initialize the computations
    %
    % Args:
    %     empty (none):       Generate default self variable assuming as products LS = Soot formation     
    %
    % Optional Args:
    %     LS (cell):          List of species
    %     obj (class):        Class combustion_toolbox_app (GUI)
    %     type (str):         Type 'fast' initialize from the given Databases
    %     DB_master (struct): Master database
    %     DB (struct) :       Database with custom thermodynamic polynomials functions generated from NASA's 9 polynomials fits
    % 
    % Returns:
    %     self (struct): Data of the mixture (initialization - empty), conditions, and databases 

    try
        [self, LS, FLAG_FAST] = initialize(varargin);
        self.E = Elements();
        self.S = Species();
        self.C = Constants();
        self.Misc = Miscellaneous();
        self.PD = ProblemDescription();
        self.PS = ProblemSolution();
        self.TN = TunningProperties();
        self = constructor(self, LS, FLAG_FAST);
        if ~nargin || ~isa(varargin{1,1}, 'combustion_toolbox_app') || ~isa(varargin{1,1}, 'combustion_toolbox_app_old') || ~isa(varargin{1,1}, 'combustion_toolbox_app_original') || (~strcmpi(varargin{1,1}, 'fast') && nargin < 4) 
            self = Initialize(self);
        end
    catch ME
      errorMessage = sprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
      ME.stack(1).name, ME.stack(1).line, ME.message);
      fprintf('%s\n', errorMessage);
      uiwait(warndlg(errorMessage));
    end
end

function [self, LS, FLAG_FAST] = initialize(varargin)
    varargin = varargin{1,1};
    nargin = length(varargin); % If varargin is empty, by default nargin will return 1, not 0.
    self = struct();
    LS = [];
    FLAG_FAST = false;
    if nargin
        if strcmpi(varargin{1,1}, 'fast')
            FLAG_FAST = true;
            self.DB_master = varargin{1,2};
            self.DB = varargin{1,3};
            if nargin == 4
                LS = varargin{1,4};
            end
            return
        end
        if isa(varargin{1,1}, 'combustion_toolbox_app') || isa(varargin{1,1}, 'combustion_toolbox_app_old') || isa(varargin{1,1}, 'combustion_toolbox_app_original')
            self = varargin{1,1};
            if nargin == 2
                LS = varargin{1,2};
            end
        else
            LS = varargin{1,1};
        end
    end
end

function self = constructor(self, LS, FLAG_FAST)
    % FLAG_GUI
    self = check_GUI(self);
    % Set Database
    FLAG_REDUCED_DB = false;
    self = set_DB(self, FLAG_REDUCED_DB, FLAG_FAST);
    % Set ListSpecies
    if ~isempty(LS)
        self = ListSpecies(self, LS);
    else
        self = ListSpecies(self);
    end
    % Set Contained elements
    self = ContainedElements(self);
    % Timer
    self.Misc.timer_0 = tic;
end

function self = check_GUI(self)
    if isa(self, 'combustion_toolbox_app') || isa(self, 'combustion_toolbox_app_old') || isa(self, 'combustion_toolbox_app_original')
        self.Misc.FLAG_GUI = true;
    end
end