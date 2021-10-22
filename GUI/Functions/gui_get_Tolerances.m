function self = gui_get_tolerances(self)
    % Get tolerance from GUI and update values
    self.C.tolN = self.TraceoptionEditField.Value;                % Tolerance of the gibbs minimization method
    self.C.tol0 = self.RootFindingMethodEditField.Value;          % Tolerance of the root finding algorithm
    self.C.itMax = self.MaxiterationsRFMEditField.Value;          % Max number of iterations - root finding method
    self.C.root_T0_l = self.RFMT0_LEditField.Value;               % First guess T[K] left branch - root finding method
    self.C.root_T0_r = self.RFMT0_REditField.Value;               % First guess T[K] right branch - root finding method
    self.C.root_T0 = self.RFMT0EditField.Value;                   % Guess T[K] if it's of previous range - root finding method
    self.C.tol_shocks = self.ShocksandDetonationsEditField.Value; % Tolerance of shocks routines
    self.C.it_shocks = self.MaxiterationsSDEditField.Value;       % Max number of iterations - shocks and detonations
end