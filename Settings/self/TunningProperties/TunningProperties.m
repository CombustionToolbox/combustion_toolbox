function self = TunningProperties()
    self.description = "Tunning properties";
    self.tolN = 1e-14;          % Tolerance of the gibbs minimization method
    self.tol_pi_e = 1e-4;       % Tolerance of the Lagrangian multiplier for ions divided by RT
    self.tol0 = 1e-3;           % Tolerance of the root finding algorithm
    self.root_method = @steff;  % Method for root finding
    self.itMax = 30;            % Max number of iterations - root finding method
    self.root_T0_l = 325;       % First guess T[K] left branch - root finding method
    self.root_T0_r = 1500;      % First guess T[K] right branch - root finding method
    self.root_T0   = 2000;      % Guess T[K] if it's of previous range - root finding method
    self.tol_shocks = 5e-5;     % Tolerance of shocks routines
    self.volumeBoundRation = 5; % Initial guess ratio shocks
    % Deprecated
    self.factor_c = 1.0;              % Tunning factor to adjust the theoretical soot equivalence ratio
    self.guess = [2000,2000,0,1.5,2]; % Initial guess shocks (deprecated - SD routines)
    self.ERRFT = 1e-6;                % Tolerances [shocks - deprecated]
    self.ERRFV = 1e-6;                % Tolerances [shocks - deprecated]
    self.ERRFU = 1e-4;                % Tolerances [detonations - deprecated]
end