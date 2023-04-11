function results = therm_effects_multicomponent(varargin)
    % Routine to perform a Linear Interaction Analysis (LIA) of
    % the interaction between a weakly turbulent free stream and a
    % hypersonic shock wave in multi-component mixtures. By default
    % considers an air mixture at sea level (T1 = 300 K and p1 = 1 atm),
    % and a wide range of pre-shock Mach numbers M1 = [1, 30].
    %
    % These calculations are based on our previous theoretical work [1]
    % and have been extended to multi-component mixtures [2] using the
    % Combustion Toolbox [3].
    %
    % Note:
    %     LIA results from [1] and [2] are computed using Mathematica
    %     and there can be some numerical differences with the results
    %     obtained with this MATLAB routine.
    %
    % References:
    %     [1] Huete, C., Cuadra, A., Vera, M., Urzay, & J. (2021). Thermochemical
    %         effects on hypersonic shock waves interacting with weak turbulence.
    %         Physics of Fluids 33, 086111 (featured article). DOI: 10.1063/5.0059948.
    %
    %     [2] Cuadra, A., Vera, M., Di Renzo, M., & Huete, C. (2023). Linear Theory
    %         of Hypersonic Shocks Interacting with Turbulence in Air. In 2023 AIAA
    %         SciTech Forum, National Harbor, USA. DOI: 10.2514/6.2023-0075.
    %
    %     [3] Cuadra, A., Huete, C., Vera, M., (2022). Combustion Toolbox:
    %         A MATLAB-GUI based open-source tool for solving gaseous
    %         combustion problems. Zenodo. DOI: 10.5281/zenodo.5554911.
    %
    % Optional Name-Value Pairs Args:
    %      * T (float): Temperature of the reactants [K]
    %      * p (float): Pressure of the reactants [bar]
    %      * species (cell): Reactants
    %      * moles (float): Molar composition of the reactants
    %      * LS (cell): List of possible products
    %      * z (float): Altitude above sea level [m]
    %      * FLAG_RESULTS (bool): Flag to show results in the command window
    %      * Npoints (int): Approximate number of points for the calculations
    %      * MAX_MACH (float): Maximum pre-shock Mach number for the calculations
    %
    % Returns:
    %      results (struct): Results of the Linear Interaction Analysis
    %
    % Examples:
    %     * therm_effects_multicomponent()
    %     * therm_effects_multicomponent('z', 0)
    %     * therm_effects_multicomponent('T', 300, 'p', 1.01325)
    
    % Default values
    filename = 'thermo_num.mat';
    FLAG_EXPORT_JUMP = false;
    FLAG_LIA = false;
    
    % Get RH jump conditions
    [R, P, T, M1, M2, Gammas] = compute_jump_conditions(varargin{:});
    
    % Export results
    if FLAG_EXPORT_JUMP
        clearvars varargin
        save(filename);
    end

    % Compute LIA
    if FLAG_LIA
        results = compute_LIA(R, M2, Gammas);
    end

end