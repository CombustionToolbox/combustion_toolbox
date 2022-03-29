function self = ProblemDescription()
    % Initialize struct with problem description data
    % 
    % Returns:
    %     self (struct): struct with problem description data

    % Description
    self.description = "Problem description";
    % Variables
    self.CompleteOrIncomplete = "incomplete";   % Default combustion assumption (deprecated)
    self.ProblemType = 'TP';                    % Default problem type
    %   * Properties matrices
    self.R_Fuel = [];                           % Fuel property matrix
    self.R_Oxidizer = [];                       % Oxidizer property matrix
    self.R_Inert = [];                          % Inert property matrix
    %   * Equivalence ratio (struct)
    self.phi.description = "Equivalence ratio"; % Description
    self.phi.value = 1.0;                       % Default value (stochiometric)
    self.phi.t = [];                            % Theoretical value for number of moles diatomic oxygen: phi = phi.t / phi.st
    %   * Fuel (struct)
    self.Fuel.x = [];                           % C atoms in the fuel mixture
    self.Fuel.y = [];                           % H atoms in the fuel mixture
    self.Fuel.z = [];                           % O atoms in the fuel mixture
    self.Fuel.w = [];                           % N atoms in the fuel mixture
    %   * Properties
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
    %   * Mixture conditions
    self.S_Fuel = [];               % Cell with the list of fuel species in the mixture
    self.N_Fuel = [];               % Vector with the number of moles of the fuel species in the mixture
    self.T_Fuel = [];               % Vector with the temperature values of the fuel species in the mixture
    self.S_Oxidizer = [];           % Cell with the list of oxidizer species in the mixture
    self.N_Oxidizer = [];           % Vector with the number of moles of the oxidizer species in the mixture
    self.T_Oxidizer = [];           % Vector with the temperature values of the oxidizer species in the mixture
    self.S_Inert = [];              % Cell with the list of inert species in the mixture
    self.N_Inert = [];              % Vector with the number of moles of the inert species in the mixture
    self.T_Inert = [];              % Vector with the temperature values of the inert species in the mixture
    self.proportion_inerts_O2 = []; % Proportion Inerts / O2 [-]
    self.ionization = false; % Flag ionization
end