function self = App(varargin)
    % Generate self variable with all the data required to initialize the computations
    %
    % Optional Args:
    %     * LS (cell): List of species
    %     * obj (class): Class combustion_toolbox_app (GUI)
    %     * type (str): If value is fast initialize from the given Databases
    %     * DB_master (struct): Master database
    %     * DB (struct): Database with custom thermodynamic polynomials functions generated from NASAs 9 polynomials fits
    %
    % Returns:
    %     self (struct): Data of the mixture (initialization - empty), conditions, and databases

    try
        [self, LS, FLAG_FAST, FLAG_COPY] = unpack(varargin);

        if ~FLAG_COPY
            self.E = Elements();
            self.S = Species();
            self.C = Constants();
            self.Misc = Miscellaneous();
            self.PD = ProblemDescription();
            self.PS = ProblemSolution();
            self.TN = TuningProperties();
        else
            self = reset_copy(self);
        end

        self = constructor(self, LS, FLAG_FAST);

        if ~nargin || ~isobject(varargin{1, 1}) || (~strcmpi(varargin{1, 1}, 'fast') && nargin < 4) || (~strcmpi(varargin{1, 1}, 'copy') && nargin < 3)
            self = initialize(self);
        end

    catch ME
        print_error(ME);
    end

end

% SUB-PASS FUNCTIONS
function [self, LS, FLAG_FAST, FLAG_COPY] = unpack(varargin)
    varargin = varargin{1, 1};
    nargin = length(varargin); % If varargin is empty, by default nargin will return 1, not 0.
    self = struct();
    LS = [];
    FLAG_FAST = false;
    FLAG_COPY = false;

    if nargin

        if strcmpi(varargin{1, 1}, 'fast')
            FLAG_FAST = true;
            self.DB_master = varargin{1, 2};
            self.DB = varargin{1, 3};

            if nargin == 4
                LS = varargin{1, 4};
            end

            return
        elseif strcmpi(varargin{1, 1}, 'copy')
            FLAG_FAST = true;
            FLAG_COPY = true;
            self = varargin{1, 2};

            if nargin == 3
                LS = varargin{1, 3};
                self.S.FLAG_COMPLETE = false;
            end

            return
        end

        if isobject(varargin{1, 1})
            self = varargin{1, 1};

            if nargin == 2
                LS = varargin{1, 2};
            end

        else
            LS = varargin{1, 1};
        end

    end

end

function self = constructor(self, LS, FLAG_FAST)
    % Timer
    self.Misc.timer_0 = tic;
    % FLAG_GUI
    self = check_GUI(self);
    % Set Database
    FLAG_REDUCED_DB = false;
    self = set_DB(self, FLAG_REDUCED_DB, FLAG_FAST);
    % Assign additional inputs (critical temperature, critical pressure, acentric factor)
    self.DB = assign_DB_PR_eos(self.DB);
    % Set ListSpecies
    if ~isempty(LS)
        self = list_species(self, LS);
    else
        self = list_species(self);
    end

    % Set Contained elements
    self = contained_elements(self);
end

function self = check_GUI(self)

    if isa(self, 'combustion_toolbox_app')
        self.Misc.FLAG_GUI = true;
    end

end

function self = reset_copy(self)
    % In case we want to solve the same base problem from a given self
    % variable, we have to reset several values

    % Reset check reactants (FOI): Fuel, Oxidizer, Inert
    self.Misc.FLAG_FOI = true;
    % Reset moles of oxidizer and inerts in case they are computed from the
    % equivalence ratio
    if ~isempty(self.PD.phi.value)
        self.PD.N_Oxidizer = [];

        if ~isempty(self.PD.ratio_inerts_O2)
            self.PD.N_Inert = [];
        end

    end

end
