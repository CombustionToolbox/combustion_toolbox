function self = TuningProperties()
    % Initialize struct with tunning properties attributes
    % 
    % Attributes:
    %     FLAG_FAST (bool):       Flag indicating use guess composition of the previous computation (default: false)
    %     itMax_gibbs (float):    Max number of iterations - Gibbs/Helmholtz minimization method    (default: 70)
    %     itMax_ions (float):     Max number of iterations - charge balance (ions)                  (default: 30)
    %     tolN (float):           Tolerance of the Gibbs/Helmholtz minimization method              (default: 1e-14)
    %     tolE (float):           Tolerance of the mass balance                                     (default: 1e-06)
    %     tol_pi_e (float):       Tolerance of the dimensionless Lagrangian multiplier - ions       (default: 1e-04)
    %     tol0 (float):           Tolerance of the root finding algorithm                           (default: 1e-03)
    %     root_method (function): Method for root finding                                           (default: newton)
    %     itMax (float):          Max number of iterations - root finding method - HP, EV, SP, SV   (default: 30)
    %     root_T0_l (float):      First guess T [K] left branch - root finding method               (default: 1000)
    %     root_T0_r (float):      First guess T [K] right branch - root finding method              (default: 3000)
    %     root_T0 (float):        Guess T[K] if it's of previous range - root finding method        (default: 3000)
    %     tol_shocks (float):     Tolerance of shocks routines                                      (default: 5e-05)
    %     it_shocks (float):      Max number of iterations - shocks and detonations                 (default: 50)
    %     Mach_thermo (float):    Preshock Mach number above which T2_guess will be computed considering h2 = h1 + u1^2 / 2 (default: 2)
    %     tol_oblique (float):    Tolerance oblique shocks                                          (default: 1e-03)
    %     it_oblique (float):     Max number of iterations - oblique shocks                         (default: 20)
    %     N_points_polar (float): Number of points to compute shock polar                           (default: 100)
    %     it_guess_det (float);   Max number of iterations - guess detonation                       (default: 5)
    %     tol_rocket (float):     Tolerance rocket performance                                      (default: 1e-04)
    %     it_rocket (float):      Max number of iterations - rocket performance                     (default: 10)
    %
    % Returns:
    %     self (struct): struct with tunning properties data


    % Description
    self.description = "Tuning properties";
    % Attributes
    % * Flags
    self.FLAG_FAST = true;     % Flag indicating use guess composition of the previous computation
    self.FLAG_TCHEM_FROZEN = false; % Flag to considerate a thermochemically frozen gas
    % * Chemical equilibrium TP, TV
    self.itMax_gibbs = 70;      % Max number of iterations - Gibbs/Helmholtz minimization method
    self.itMax_ions = 30;       % Max number of iterations - charge balance (ions)
    self.tolN = 1e-14;          % Tolerance of the composition of the mixture
    self.tolN_guess = 1e-6;     % Tolerance of the composition of the mixture (guess)
    self.tol_gibbs = 1e-5;      % Tolerance of the Gibbs/Helmholtz minimization method
    self.tolE = 1e-6;           % Tolerance of the mass balance
    self.tol_pi_e = 1e-4;       % Tolerance of the dimensionless Lagrangian multiplier - ions
    self.T_ions = 500;          % Minimum temperature [K] to consider ionized species
    % * Chemical equilibrium HP, EV, SP, SV
    self.tol0 = 1e-3;           % Tolerance of the root finding algorithm
    self.root_method = @newton; % Root finding method [newton (2nd order), steff (2nd order), or nsteff (3rd order)]
    self.itMax = 30;            % Max number of iterations - root finding method
    self.root_T0_l = 1000;      % First guess T[K] left branch - root finding method
    self.root_T0_r = 3000;      % First guess T[K] right branch - root finding method
    self.root_T0   = 3000;      % Guess T[K] if it's of previous range - root finding method
    % * Shocks and detonations 
    self.tol_shocks = 1e-5;     % Tolerance of shocks/detonations kernel
    self.it_shocks = 50;        % Max number of iterations - shocks and detonations
    self.Mach_thermo = 2;       % Preshock Mach number above which T2_guess will be computed considering h2 = h1 + u1^2 / 2
    self.tol_oblique = 1e-3;    % Tolerance oblique shocks algorithm
    self.it_oblique = 20;       % Max number of iterations - oblique shocks
    self.N_points_polar = 100;  % Number of points to compute shock/detonation polar curves
    self.tol_limitRR = 1e-4;    % Tolerance to calculate the limit of regular reflections
    self.it_limitRR = 10;       % Max number of iterations - limit of regular reflections
    self.it_guess_det = 5;      % Max number of iterations - guess detonation
    % * Rocket propellant performance
    self.tol_rocket = 1e-4;     % Tolerance rocket performance
    self.it_rocket = 10;        % Max number of iterations - rocket performance
    % * Equation of State (EoS)
    self.tol_eos = 1e-4;        % Tolerance equation of state (EoS) - only Newton-Raphson method (deprecated)
    self.it_eos = 30;           % Max number of iterations - equation of state (EoS) - only Newton-Raphson method (deprecated)
end