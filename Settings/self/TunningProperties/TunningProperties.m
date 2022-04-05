function self = TunningProperties()
    % Initialize struct with tunning properties data
    % 
    % Returns:
    %     self (struct): struct with tunning properties data

    % Description
    self.description = "Tunning properties";
    % Variables
    %   * Flags
    self.FLAG_FAST = false;     % Flag indicating use guess composition of the previous computation
    %   * Chemical equilibrium TP, TV
    self.tolN = 1e-15;          % Tolerance of the gibbs minimization method
    self.tolE = 1e-6;           % Tolerance of the mass balance
    self.tol_pi_e = 1e-4;       % Tolerance of the Lagrangian multiplier for ions divided by RT
    %   * Chemical equilibrium HP, EV, SP, SV
    self.tol0 = 1e-3;           % Tolerance of the root finding algorithm
    self.root_method = @newton; % Method for root finding
    self.itMax = 30;            % Max number of iterations - root finding method
    self.root_T0_l = 300;       % First guess T[K] left branch - root finding method
    self.root_T0_r = 1500;      % First guess T[K] right branch - root finding method
    self.root_T0   = 2000;      % Guess T[K] if it's of previous range - root finding method
    %   * Shocks and detonations 
    self.tol_oblique = 1e-3;    % Tolerance oblique shocks
    self.it_oblique = 20;       % Max number of iterations - oblique shocks
    self.N_points_polar = 100;  % Number of points to compute shock polar
    self.tol_shocks = 5e-5;     % Tolerance of shocks routines
    self.it_shocks = 50;        % Max number of iterations - shocks and detonations
    self.volumeBoundRation = 5; % Initial guess ratio shocks
    %   * Rocket propellant performance
    self.tol_rocket = 0.4*1e-4; % Tolerance rocket performance
    self.it_rocket = 2;         % Max number of iterations - rocket performance
end