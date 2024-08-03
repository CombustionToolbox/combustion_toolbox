function self = run_CT(varargin)
    % A generalized function to run Combustion Toolbox for a given set of
    % inputs. Otherwise, it will run the predefined case.
    
    % cd('../../');

    % DEFAULT VALUES
    DB = [];
    DB_master = [];
    species = [];
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
    problemType = 'HP';
    tolN = 1e-18;
    N_points_polar = 100;
    FLAG_PHI = true;
    FLAG_IAC = true;
    FLAG_TCHEM_FROZEN = false;
    FLAG_FROZEN = false;
    FLAG_FAST = true;
    FLAG_RESULTS = false;
    vP_vR = [];
    Aratio_c = [];
    Aratio = [];
    % GET INPUTS
    for i = 1:2:nargin

        switch lower(varargin{i})
            case {'db'}
                DB = varargin{i + 1};
            case {'problemtype', 'problem'}
                problemType = varargin{i + 1};
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
            case {'vp_vr'}
                vP_vR = varargin{i + 1};
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
            case {'flag_tchem_frozen', 'flag_tchem'}
                FLAG_TCHEM_FROZEN = varargin{i + 1};
            case 'flag_frozen'
                FLAG_FROZEN = varargin{i + 1};
            case 'aratio_c'
                Aratio_c = varargin{i + 1};
            case 'aratio'
                Aratio = varargin{i + 1};
            case 'flag_fast'
                FLAG_FAST = varargin{i + 1};
            case 'flag_results'
                FLAG_RESULTS = varargin{i + 1};
        end

    end

    % Initialiaze
    if isempty(DB)
        DB = combustiontoolbox.databases.NasaDatabase();
    end
    
    system = combustiontoolbox.core.ChemicalSystem(DB, species);

    % MISCELLANEOUS
    self.Misc.FLAG_RESULTS = FLAG_RESULTS;

    % Define the initial mixture composition
    mix1 = combustiontoolbox.core.Mixture(system);
    set(mix1, S_Fuel, N_Fuel);
    % set(mix1, S_Oxidizer, 'oxidizer', N_Oxidizer);
    % set(mix1, S_Inert, 'inert', N_Inert);
    
    % Define the thermodynamic state
    mix1Array = setProperties(mix1, 'temperature', Temp_mix1 * ones(size(Temp_mix2)), 'pressure', Pressure_mix1);
    
    % Select solver
    switch upper(problemType)
        case {'TP', 'HP', 'SP', 'TV', 'EV', 'SV'}
            solver = combustiontoolbox.equilibrium.EquilibriumSolver('problemType', problemType);
    end
    
    % Solve chemical transformation at equilibrium
    % tic
    numCases = length(Temp_mix2);
    mix2Array = setProperties(mix1, 'temperature', Temp_mix2, 'pressure', Pressure_mix2);
    solver.solveArray(mix2Array);


    % % TUNNING PROPERTIES
    % self.TN.tolN = tolN;
    % self.TN.N_points_polar = N_points_polar;
    % self.TN.FLAG_FAST = FLAG_FAST;
    % 
    % % INITIAL CONDITIONS
    % self = set_prop(self, 'TR', Temp_mix1, 'pR', Pressure_mix1);
    % 
    % if FLAG_PHI
    %     self = set_prop(self, 'phi', EquivalenceRatio);
    % end
    % 
    % if ~isempty(vP_vR)
    %     self = set_prop(self, 'vP_vR', vP_vR);
    % end
    % 
    % if ~isempty(Aratio_c)
    %     self = set_prop(self, 'Aratio_c', Aratio_c);
    % end
    % 
    % if ~isempty(Aratio)
    %     self = set_prop(self, 'Aratio', Aratio);
    % end
    % 
    % self.PD.FLAG_IAC = FLAG_IAC;
    % self.PD.FLAG_TCHEM_FROZEN = FLAG_TCHEM_FROZEN;
    % self.PD.FLAG_FROZEN = FLAG_FROZEN;
    % self.PD.S_Fuel = S_Fuel; self.PD.N_Fuel = N_Fuel;
    % self.PD.S_Oxidizer = S_Oxidizer; self.PD.N_Oxidizer = N_Oxidizer;
    % self.PD.S_Inert = S_Inert; self.PD.N_Inert = N_Inert;
    % self.PD.ratio_oxidizers_O2 = ratio_oxidizers_O2;
    % self.PD.ratio_inerts_O2 = ratio_inerts_O2;
    % % ADDITIONAL INPUTS (DEPENDS OF THE PROBLEM SELECTED)
    % self = set_prop(self, 'pP', Pressure_mix2, 'TP', Temp_mix2);
    % 
    % if exist('Velocity', 'var')
    %     self = set_prop(self, 'u1', Velocity, 'phi', 1 * ones(1, length(Velocity)));
    % end
    % 
    % % SOLVE PROBLEM
    % self = solve_problem(self, problemType);

    % SET PROPERTIES AS CEA
    for i = numCases:-1:1
        % Mix 1
        self.PS.strR{i}.h = enthalpy_mass(mix1Array(i));
        self.PS.strR{i}.e = intEnergy_mass(mix1Array(i));
        self.PS.strR{i}.g = gibbs_mass(mix1Array(i));
        self.PS.strR{i}.S = entropy_mass(mix1Array(i));
        self.PS.strR{i}.cP = cp_mass(mix1Array(i));
        self.PS.strR{i}.cV = mix1Array(i).cp / mix1Array(i).gamma_s;
        % Mix 2
        self.PS.strP{i}.h = enthalpy_mass(mix2Array(i));
        self.PS.strP{i}.e = intEnergy_mass(mix2Array(i));
        self.PS.strP{i}.g = gibbs_mass(mix2Array(i));
        self.PS.strP{i}.S = entropy_mass(mix2Array(i));
        self.PS.strP{i}.cP = cp_mass(mix2Array(i));
        self.PS.strP{i}.cV = mix2Array(i).cp / mix2Array(i).gamma_s;

        % if contains(problemType, '_R')
        %     self.PS.strP{i}.u = self.PS.strR{i}.u;
        % end
        % 
        % if contains(problemType, 'ROCKET')
        %     self.PS.mix2_c{i}.h = enthalpy_mass(self.PS.mix2_c{i});
        %     self.PS.mix2_c{i}.e = intEnergy_mass(self.PS.mix2_c{i});
        %     self.PS.mix2_c{i}.g = gibbs_mass(self.PS.mix2_c{i});
        %     self.PS.mix2_c{i}.S = entropy_mass(self.PS.mix2_c{i});
        %     self.PS.mix2_c{i}.cP = cp_mass(self.PS.mix2_c{i});
        %     self.PS.mix2_c{i}.cV = self.PS.mix2_c{i}.cP / self.PS.mix2_c{i}.gamma_s;
        % 
        %     self.PS.mix3{i}.h = enthalpy_mass(self.PS.mix3{i});
        %     self.PS.mix3{i}.e = intEnergy_mass(self.PS.mix3{i});
        %     self.PS.mix3{i}.g = gibbs_mass(self.PS.mix3{i});
        %     self.PS.mix3{i}.S = entropy_mass(self.PS.mix3{i});
        %     self.PS.mix3{i}.cP = cp_mass(self.PS.mix3{i});
        %     self.PS.mix3{i}.cV = self.PS.mix3{i}.cP / self.PS.mix3{i}.gamma_s;
        % end
        % 
        % if contains(problemType, 'SHOCK') || contains(problemType, 'DET')
        %     self.PS.strR{i}.u_preshock = self.PS.strR{i}.u;
        %     self.PS.strR{i}.u_postshock = self.PS.strP{i}.v_shock;
        %     self.PS.strP{i}.u_preshock = self.PS.strR{i}.u;
        %     self.PS.strP{i}.u_postshock = self.PS.strP{i}.v_shock;
        % end

    end

end
