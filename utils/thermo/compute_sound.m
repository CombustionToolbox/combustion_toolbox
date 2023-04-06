function sound = compute_sound(T, p, species, composition, varargin)
    % Routine to compute sound velocity [m/s] for a given
    % temperature-pressure profile
    %
    % Args:
    %     T (float): Temperature [K]
    %     p (float): Pressure [bar]
    %     species (cell): List of species
    %     composition (float): composition of species (mol) 
    %
    % Optional Name-Value Pairs Args:
    %     self (struct): Data of the mixtures, conditions, databases
    %
    % Returns:
    %     sound (float): Sound velocity [m/s]
    %
    % Examples:
    %     sound = compute_sound(300, 1, {'H2', 'O2'}, [1, 1])
    %     sound = compute_sound(300, 1, {'H2', 'O2'}, [1, 1], 'self', self)
    
    % Default
    FLAG_SELF = false;

    % Unpack
    for i = 1:2:nargin - 4
        switch lower(varargin{i})
            case 'self'
                self = varargin{i + 1};
                FLAG_SELF = true;
        end
    end
    
    % Check inputs
    if ~iscell(species)
        species = {species};
    end

    % Initialization
    if FLAG_SELF
        self = App('FAST',  self.DB_master, self.DB, species);
    else
        self = App(species);
    end

    % Set mixture
    self.PD.S_Fuel = species;
    self.PD.N_Fuel = composition / sum(composition);
    
    % Set temperature-pressure profile
    self = set_prop(self, 'TR', T, 'pR', p);
    
    % Set Problem Type
    self.PD.ProblemType = 'None';
    
    % Check inputs and set length of the loop
    self = check_inputs(self);

    for i = self.C.l_phi:-1:1
        % Set problem conditions by case
        self = set_problem_conditions(self, i);
        % Compute properties matrix
        self.PD.R_Fuel = set_species(self, self.PD.S_Fuel, self.PD.N_Fuel, self.PD.TR.value);
        % Compute properties mixture
        self.PS.strR_Fuel = compute_properties(self, self.PD.R_Fuel, self.PD.pR.value, self.PD.TR.value);
        % Get sound velocity
        sound(i) = soundspeed(self.PS.strR_Fuel);
    end

end

% SUB-PASS FUNCTIONS
function self = set_problem_conditions(self, i)
    % Set problem conditions per case
    if ~isfield(self.PD, 'range_name')
        return
    end

    if ~strcmpi(self.PD.range_name, 'phi')
        self.PD.(self.PD.range_name).value = self.PD.range(i);
    end

end