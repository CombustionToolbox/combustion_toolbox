function self = ProblemDescription()
    % Initialize struct with problem description data
    % 
    % Returns:
    %     self (struct): struct with problem description data

    % Description
    self.description = "Problem description";
    % Attributes:
    self.CompleteOrIncomplete = "incomplete";   % Default combustion assumption (deprecated)
    self.ProblemType = 'TP';                    % Default problem type
    % * Properties matrices
    self.R_Fuel = [];                           % Fuel property matrix
    self.R_Oxidizer = [];                       % Oxidizer property matrix
    self.R_Inert = [];                          % Inert property matrix
    % * Equivalence ratio (struct)
    self.phi.description = "Equivalence ratio"; % Description
    self.phi.value = 1.0;                       % Default value (stochiometric)
    self.phi.t = [];                            % Theoretical value for number of moles diatomic oxygen: phi = phi.t / phi.st
    % * Fuel (struct)
    self.Fuel.x = [];                           % C atoms in the fuel mixture
    self.Fuel.y = [];                           % H atoms in the fuel mixture
    self.Fuel.z = [];                           % O atoms in the fuel mixture
    self.Fuel.w = [];                           % N atoms in the fuel mixture
    % * Properties
    self.TR.description = "Temperature of reactants";
    self.TR.value = []; % [K]
    self.pR.description = "Pressure of reactants";
    self.pR.value = []; % [bar]
    self.TP.description = "Temperature of products";
    self.TP.value = []; % [K]
    self.pP.description = "Pressure of products";
    self.pP.value = []; % [bar]
    self.vP_vR.description = "Volume relation Products/Reactants";
    self.vP_vR.value = []; % [-]
    self.u1.description = "Incident shock velocity";
    self.u1.value = []; % [m/s]
    self.overdriven.description = "Overdriven shock velocity";
    self.overdriven.value = []; % [-]
    self.theta.description = "Deflection angle - oblique shocks";
    self.theta.value = []; % [deg]
    self.beta.description = "Wave angle - oblique shocks";
    self.beta.value = []; % [deg]
    self.Aratio.description = "Area ratio exit/throat - rocket";
    self.Aratio.value = []; % [-]
    self.Aratio_c.description = "Area ratio combustion chamber/thoat - rocket";
    self.Aratio_c.value = []; % [-]
    % * Mixture conditions
    self.S_Fuel = [];             % Cell with the list of fuel species in the mixture
    self.N_Fuel = [];             % Vector with the number of moles of the fuel species in the mixture
    self.T_Fuel = [];             % Vector with the temperature values of the fuel species in the mixture
    self.S_Oxidizer = [];         % Cell with the list of oxidizer species in the mixture
    self.N_Oxidizer = [];         % Vector with the number of moles of the oxidizer species in the mixture
    self.T_Oxidizer = [];         % Vector with the temperature values of the oxidizer species in the mixture
    self.S_Inert = [];            % Cell with the list of inert species in the mixture
    self.N_Inert = [];            % Vector with the number of moles of the inert species in the mixture
    self.T_Inert = [];            % Vector with the temperature values of the inert species in the mixture
    self.ratio_oxidizers_O2 = 1;  % Ratio oxidizers / O2 [% moles]
    self.ratio_inerts_O2 = [];    % Ratio oxidizers / Inerts [% moles]
    self.wt_ratio_oxidizers = []; % Weight ratio percentage of oxidizer species
    self.wt_ratio_inerts = [];    % Weight ratio percentage of inert species
    % * Equation of States
    self.EOS.pressure = @eos_ideal_p; % Equation of State to compute pressure [Pa]
    self.EOS.volume = @eos_ideal; % Equation of State to compute molar volume [m3/mol]
    self.EOS.chemical_potential_imp = @mu_imp_ideal; % Compute non ideal contribution of the chemical potential (depends of the Equation of State) [J/mol]
    % * Flags
    self.FLAG_ION = false;       % Flag ionized species in the system
    self.FLAG_IAC = true;        % Flag use IAC model for rocket computations
    self.FLAG_SUBSONIC = false;  % Flag subsonic Area ratio
    self.FLAG_EOS = false;       % Flag to use non-ideal Equation of States (EoS)
end