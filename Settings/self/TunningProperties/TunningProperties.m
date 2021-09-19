function self = TunningProperties()
    self.description = "Tunning properties";
    self.factor_c = 1.0;
    % SHOCK and DET code
    self.guess = [2000,2000,0,1.5,2]; % Initial guess
    self.ERRFT = 1e-6; % Tolerances
    self.ERRFU = 1e-4; % Tolerances
    self.ERRFV = 1e-6; % Tolerances
    self.itMax = 30;
    self.tolN = 1e-14; % Tolerance of the gibbs minimization method
    self.tol0 = 1e-3;  % Tolerance of the root finding algorithm
    self.root_method = @steff; % Method for root finding
    self.tol_shocks = 5e-5;  % Tolerance of shocks routines
    self.volumeBoundRation = 5; % Initial guess ratio shocks
end