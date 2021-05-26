function self = TunningProperties()
    self.description = "Tunning properties";
    self.factor_c = 1.0;
    % SHOCK and DET code
    self.guess = [2000,2000,0,1.5,2]; % Initial guess
    self.ERRFT = 1e-4; % Tolerances
    self.ERRFU = 1e-4; % Tolerances
    self.ERRFV = 1e-4; % Tolerances
end