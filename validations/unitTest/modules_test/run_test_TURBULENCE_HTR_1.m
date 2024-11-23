function [max_rel_error_decomposition, max_rel_error_spectra] = run_test_TURBULENCE_HTR_1()
    
    %
    % Definitions
    filename = 'velocityField.h5';
    filePath = which(filename);

    % Load data
    rho = double(h5read(filePath, ['/', 'density']));
    u = double(h5read(filePath, ['/', 'velocity/u']));
    v = double(h5read(filePath, ['/', 'velocity/v']));
    w = double(h5read(filePath, ['/', 'velocity/w']));

    % Convert velocity to VelocityField object
    velocityField = VelocityField(u, v, w);
    
    % Initialize HelmholtzSolver
    solver = HelmholtzSolver();
    
    % Perform Helmholtz decomposition
    [solenoidal, dilatational] = solver.solve(velocityField, 'rho', rho);
    
    % Analyze turbulence spectra using TurbulenceSpectra
    analyzer = TurbulenceSpectra();
    
    % Compute energy spectra for the original, solenoidal, and dilatational fields
    [EK1, k1] = analyzer.getEnergySpectra(velocityField); % Original field
    [EK2, k2] = analyzer.getEnergySpectra(solenoidal);    % Solenoidal component
    [EK3, k3] = analyzer.getEnergySpectra(dilatational);  % Dilatational component

    % Compare results
    
end