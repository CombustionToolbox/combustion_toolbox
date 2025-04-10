function [max_rel_error_decomposition, max_rel_error_spectra] = run_test_TURBULENCE_HTR_1()
    
    % Import packages
    import combustiontoolbox.turbulence.*

    % Definitions
    filename = 'velocityField.h5';
    filePath = which(filename);

    % Load data
    rho = double(h5read(filePath, ['/', 'density']));
    u = double(h5read(filePath, ['/', 'velocity/u']));
    v = double(h5read(filePath, ['/', 'velocity/v']));
    w = double(h5read(filePath, ['/', 'velocity/w']));

    % Convert velocity to VelocityField object
    velocity = VelocityField(u, v, w);
    
    % Initialize HelmholtzSolver
    solver = HelmholtzSolver();
    
    % Perform Helmholtz decomposition
    [solenoidal, dilatational, velocity] = solver.solve(velocity, 'density', rho);
    
    % Analyze turbulence spectra using TurbulenceSpectra
    analyzer = TurbulenceSpectra();
    
    % Compute energy spectra for the original, solenoidal, and dilatational fields
    [EK1, k1] = analyzer.getEnergySpectraVelocity(velocity);      % Original field
    [EK2, k2] = analyzer.getEnergySpectraVelocity(solenoidal);    % Solenoidal component
    [EK3, k3] = analyzer.getEnergySpectraVelocity(dilatational);  % Dilatational component

    % Compare results
    figure; hold on; set(gca, 'XScale', 'log');
    plot(k1, EK1(1:length(k1)), 'k-', 'LineWidth', 1.6);
    plot(k2, EK2(1:length(k2)), 'k--', 'LineWidth', 1.6);
    plot(k3, EK3(1:length(k3)), 'k:', 'LineWidth', 1.6);
end