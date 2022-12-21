function self = run_CT(varargin)
    % A generalized function to run Combustion Toolbox for a given set of
    % inputs. Otherwise, it will run the predefined case.

    % DEFAULT VALUES
    species = 'Soot Formation';
    Temp_mix1 = 300;
    Temp_mix2 = 2500;
    Pressure_mix1 = 1;
    Pressure_mix2 = 1;
    EquivalenceRatio = 0.5:0.01:4;
    S_Fuel = []; N_Fuel = [];
    S_Oxidizer = {'N2', 'O2'}; N_Oxidizer = [];
    S_Inert = []; N_Inert = [];
    ratio_oxidizers_O2 = [79, 21] / 21;
    ratio_inerts_O2 = [];
    ProblemType = 'HP';
    tolN = 1e-18;
    N_points_polar = 100;
    FLAG_PHI = true;
    FLAG_IAC = true;
    Aratio_c = [];
    Aratio = [];
    % GET INPUTS
    for i = 1:2:nargin

        switch lower(varargin{i})
            case {'problemtype', 'problem'}
                ProblemType = varargin{i + 1};
            case {'listspecies', 'species'}
                species = varargin{i + 1};
            case {'tr'}
                Temp_mix1 = varargin{i + 1};
            case {'temperature', 'temp', 'tp'}
                Temp_mix2 = varargin{i + 1};
            case {'pr'}
                Pressure_mix1 = varargin{i + 1};
            case {'pressure', 'pp'}
                Pressure_mix2 = varargin{i + 1};
            case {'equivalenceratio', 'phi'}
                EquivalenceRatio = varargin{i + 1};
            case {'velocity', 'u1'}
                Velocity = varargin{i + 1};
            case {'fuel', 's_fuel'}
                S_Fuel = varargin{i + 1};

                if ~iscell(S_Fuel) && ~isempty(S_Fuel)
                    S_Fuel = {S_Fuel};
                end

            case {'oxidizer', 's_oxidizer'}
                S_Oxidizer = varargin{i + 1};

                if ~iscell(S_Oxidizer) && ~isempty(S_Oxidizer)
                    S_Oxidizer = {S_Oxidizer};
                end

            case {'inert', 's_inert'}
                S_Inert = varargin{i + 1};

                if ~iscell(S_Inert) && ~isempty(S_Inert)
                    S_Inert = {S_Inert};
                end

            case 'n_fuel'
                N_Fuel = varargin{i + 1};
                FLAG_PHI = false;
            case 'n_oxidizer'
                N_Oxidizer = varargin{i + 1};
                FLAG_PHI = false;
            case 'n_inert'
                N_Inert = varargin{i + 1};
                FLAG_PHI = false;
            case 'ratio_oxidizers_o2'
                ratio_oxidizers_O2 = varargin{i + 1};
            case 'ratio_inerts_o2'
                ratio_inerts_O2 = varargin{i + 1};
            case 'toln'
                tolN = varargin{i + 1};
            case 'n_points_polar'
                N_points_polar = varargin{i + 1};
            case 'flag_iac'
                FLAG_IAC = varargin{i + 1};
            case 'aratio_c'
                Aratio_c = varargin{i + 1};
            case 'aratio'
                Aratio = varargin{i + 1};
        end

    end

    % INITIALIZE
    self = App(species);
    % MISCELLANEOUS
    self.Misc.FLAG_RESULTS = false;
    % TUNNING PROPERTIES
    self.TN.tolN = tolN;
    self.TN.N_points_polar = N_points_polar;
    % INITIAL CONDITIONS
    self = set_prop(self, 'TR', Temp_mix1, 'pR', Pressure_mix1);

    if FLAG_PHI
        self = set_prop(self, 'phi', EquivalenceRatio);
    end

    if ~isempty(Aratio_c)
        self = set_prop(self, 'Aratio_c', Aratio_c);
    end

    if ~isempty(Aratio)
        self = set_prop(self, 'Aratio', Aratio);
    end

    self.PD.FLAG_IAC = FLAG_IAC;
    self.PD.S_Fuel = S_Fuel; self.PD.N_Fuel = N_Fuel;
    self.PD.S_Oxidizer = S_Oxidizer; self.PD.N_Oxidizer = N_Oxidizer;
    self.PD.S_Inert = S_Inert; self.PD.N_Inert = N_Inert;
    self.PD.ratio_oxidizers_O2 = ratio_oxidizers_O2;
    self.PD.ratio_inerts_O2 = ratio_inerts_O2;
    % ADDITIONAL INPUTS (DEPENDS OF THE PROBLEM SELECTED)
    self = set_prop(self, 'pP', Pressure_mix2, 'TP', Temp_mix2);

    if exist('Velocity', 'var')
        self = set_prop(self, 'u1', Velocity, 'phi', 1 * ones(1, length(Velocity)));
    end

    % INITIALIZE TIMER
    self.Misc.timer_loop = tic;
    % SOLVE PROBLEM
    self = solve_problem(self, ProblemType);
    % SAVE TIMER (COMPUTATIONAL TIME [s])
    self.Misc.timer_loop = toc(self.Misc.timer_0);
    % SET PROPERTIES AS CEA
    for i = length(self.PD.range):-1:1
        % Mix 1
        self.PS.strR{i}.h = enthalpy_mass(self.PS.strR{i});
        self.PS.strR{i}.e = intEnergy_mass(self.PS.strR{i});
        self.PS.strR{i}.g = gibbs_mass(self.PS.strR{i});
        self.PS.strR{i}.S = entropy_mass(self.PS.strR{i});
        self.PS.strR{i}.cP = cp_mass(self.PS.strR{i});
        self.PS.strR{i}.cV = self.PS.strR{i}.cP / self.PS.strR{i}.gamma_s;
        % Mix 2
        self.PS.strP{i}.h = enthalpy_mass(self.PS.strP{i});
        self.PS.strP{i}.e = intEnergy_mass(self.PS.strP{i});
        self.PS.strP{i}.g = gibbs_mass(self.PS.strP{i});
        self.PS.strP{i}.S = entropy_mass(self.PS.strP{i});
        self.PS.strP{i}.cP = cp_mass(self.PS.strP{i});
        self.PS.strP{i}.cV = self.PS.strP{i}.cP / self.PS.strP{i}.gamma_s;

        if contains(ProblemType, '_R')
            self.PS.strP{i}.u = self.PS.strR{i}.u;
        end

        if contains(ProblemType, 'ROCKET')
            self.PS.mix2_c{i}.h = enthalpy_mass(self.PS.mix2_c{i});
            self.PS.mix2_c{i}.e = intEnergy_mass(self.PS.mix2_c{i});
            self.PS.mix2_c{i}.g = gibbs_mass(self.PS.mix2_c{i});
            self.PS.mix2_c{i}.S = entropy_mass(self.PS.mix2_c{i});
            self.PS.mix2_c{i}.cP = cp_mass(self.PS.mix2_c{i});
            self.PS.mix2_c{i}.cV = self.PS.mix2_c{i}.cP / self.PS.mix2_c{i}.gamma_s;

            self.PS.mix3{i}.h = enthalpy_mass(self.PS.mix3{i});
            self.PS.mix3{i}.e = intEnergy_mass(self.PS.mix3{i});
            self.PS.mix3{i}.g = gibbs_mass(self.PS.mix3{i});
            self.PS.mix3{i}.S = entropy_mass(self.PS.mix3{i});
            self.PS.mix3{i}.cP = cp_mass(self.PS.mix3{i});
            self.PS.mix3{i}.cV = self.PS.mix3{i}.cP / self.PS.mix3{i}.gamma_s;
        end

        if contains(ProblemType, 'SHOCK') || contains(ProblemType, 'DET')
            self.PS.strR{i}.u_preshock = self.PS.strR{i}.u;
            self.PS.strR{i}.u_postshock = self.PS.strP{i}.v_shock;
            self.PS.strP{i}.u_preshock = self.PS.strR{i}.u;
            self.PS.strP{i}.u_postshock = self.PS.strP{i}.v_shock;
        end

    end

end
