function app = gui_get_tolerances(app)
    % Get tolerance from GUI and update values
    app.C.tolN = app.TraceoptionEditField.Value;                % Tolerance of the gibbs minimization method
    app.C.tol0 = app.RootFindingMethodEditField.Value;          % Tolerance of the root finding algorithm
    app.C.itMax = app.MaxiterationsRFMEditField.Value;          % Max number of iterations - root finding method
    app.C.root_T0_l = app.RFMT0_LEditField.Value;               % First guess T[K] left branch - root finding method
    app.C.root_T0_r = app.RFMT0_REditField.Value;               % First guess T[K] right branch - root finding method
    app.C.root_T0 = app.RFMT0EditField.Value;                   % Guess T[K] if it's of previous range - root finding method
    app.C.tol_shocks = app.ShocksandDetonationsEditField.Value; % Tolerance of shocks routines
    app.C.it_shocks = app.MaxiterationsSDEditField.Value;       % Max number of iterations - shocks and detonations
end