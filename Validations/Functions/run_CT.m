function self = run_CT(varargin)
    % A generalized function to run Combustion Toolbox for a given set of
    % inputs. Otherwise, it will run the predefined case.

    % DEFAULT VALUES
    species = 'Soot Formation';
    Temp = 300;
    Pressure = 1.01325;
    EquivalenceRatio = 0.5:0.01:5;
    S_Fuel = {'CH4'};
    S_Oxidizer = {'O2'};
    S_Inert    = {'N2'};
    proportion_inerts_O2 = 79/21;
    ProblemType = 'HP';
    tolN = 1e-16;
    % GET INPUTS
    for i=1:2:nargin
        switch lower(varargin{i})
            case {'problemtype', 'problem'}
                ProblemType = varargin{i+1};
            case {'listspecies', 'species'}
                species = varargin{i+1};
            case {'temperature', 'temp', 'tr'}
                Temp = varargin{i+1};
            case {'pressure', 'pr'}
                Pressure = varargin{i+1};
            case {'equivalenceratio', 'phi'}
                EquivalenceRatio = varargin{i+1};
            case {'velocity', 'u1'}
                Velocity = varargin{i+1};
            case {'fuel', 's_fuel'}
                S_Fuel = varargin{i+1};
                if ~iscell(S_Fuel) && ~isempty(S_Fuel)
                    S_Fuel = {S_Fuel};
                end
            case {'oxidizer', 's_oxidizer'}
                S_Oxidizer = varargin{i+1};
                if ~iscell(S_Oxidizer) && ~isempty(S_Oxidizer)
                    S_Oxidizer = {S_Oxidizer};
                end
            case {'inert', 's_inert'}
                S_Inert = varargin{i+1};
                if ~iscell(S_Inert) && ~isempty(S_Inert)
                    S_Inert = {S_Inert};
                end
            case 'proportion_inerts_o2'
                proportion_inerts_O2 = varargin{i+1};
            case 'toln'
                tolN = varargin{i+1};
        end
    end
    % INITIALIZE
    self = App(species);
    % MISCELLANEOUS
    self.Misc.FLAG_RESULTS = false;
    % TUNNING PROPERTIES
    self.TN.tolN = tolN;
    % INITIAL CONDITIONS
    self = set_prop(self, 'TR', Temp, 'pR', 1 * Pressure, 'phi', EquivalenceRatio);
    self.PD.S_Fuel     = S_Fuel;
    self.PD.S_Oxidizer = S_Oxidizer;
    self.PD.S_Inert    = S_Inert;
    self.PD.proportion_inerts_O2 = proportion_inerts_O2;
    % ADDITIONAL INPUTS (DEPENDS OF THE PROBLEM SELECTED)
    self = set_prop(self, 'TP', Temp, 'pP', Pressure);
    if exist('Velocity', 'var')
        self = set_prop(self, 'u1', Velocity, 'phi', 1 * ones(1, length(Velocity)));
    end
    % SOLVE PROBLEM
    self = SolveProblem(self, ProblemType);
end