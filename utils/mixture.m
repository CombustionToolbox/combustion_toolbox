function mix = mixture(T, p, species, moles, varargin)
    % Compute the properties of a mixture at a given temperature and pressure
    %
    % Args:
    %     T (float): Temperature [K]
    %     p (float): Pressure [bar]
    %     species (cell): List of species
    %     moles (cell): List of moles of each species
    %
    % Optional Name-Value Pairs Args:
    %     * phi (float): Equivalence ratio [-]
    %
    % Returns:
    %     mix (struct): Mixture properties
    %
    % Examples:
    %     * mix = mixture(300, 1, {'CH4', 'O2', 'N2'}, [1, 2, 7.52])
    %     * mix = mixture(300, 1, {'CH4', 'O2', 'N2'}, [1, 2, 7.52], 'self', self)
    %     * mix = mixture(300, 1, {'CH4', 'O2', 'N2'}, [1, 2, 7.52], 'DB_master', DB_master, 'DB', DB)

    % Defaults
    FLAG_FAST = false;

    % Unpack
    for i = 1:2:nargin-4
        switch lower(varargin{i})
            case 'self'
                self = varargin{i + 1};
                DB_master = self.DB_master;
                DB = self.DB;
                FLAG_FAST = true;
            case 'db_master'
                DB_master = varargin{i + 1};
            case 'db'
                DB = varargin{i + 1};
                FLAG_FAST = true;
        end

    end

    % Initialization
    if FLAG_FAST
        self = App('fast', DB_master, DB, species);
    else
        self = App(species);
    end

    % Set the properties of the mixture
    self = set_prop(self, 'TR', T, 'pR', p);

    % Set composition of the mixture
    self.PD.S_Oxidizer = species;
    self.PD.N_Oxidizer = moles;

    % Get Flags and length of the loop
    self = get_FLAG_N(self);

    % Compute the properties of the initial mixture
    self = define_FOI(self, 1);

    % Get mixture
    mix = self.PS.strR{1};
end