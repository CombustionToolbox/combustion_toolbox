function [R, P, T, M1, M2, Gammas] = compute_jump_conditions(varargin)
    % Compute jump conditions of non-reactive shocks. By default
    % considers an air mixture at sea level (T1 = 300 K and p1 = 1 atm),
    % and a wide range of pre-shock Mach numbers M1 = [1, 30].
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
    %      Tuple containing
    %
    %      * R (float): Density ratio rho_2 / rho_1
    %      * P (float): Pressure ratio p2 / p1
    %      * T (float): Temperature ratio T2 / T1
    %      * M1 (float): Pre-shock Mach number
    %      * M2 (float): Post-shock Mach number
    %      * Gammas (float): Inverse normalized Hugonitor slope, see Eq. (22) in [1] (Eq. (4) in [2])
    %
    % Examples:
    %      * [R, P, T, M1, M2, Gammas] = compute_jump_conditions()
    %      * [R, P, T, M1, M2, Gammas] = compute_jump_conditions('z', 0)
    %      * [R, P, T, M1, M2, Gammas] = compute_jump_conditions('T', 300, 'p', 1.01325)
    %
    % Note: 
    %      This script may use "Standard Atmosphere Functions" - FileExchange [3]
    %
    % References:
    %     [1] Huete, C., Cuadra, A., Vera, M., & Urzay, J. (2021). Thermochemical
    %         effects on hypersonic shock waves interacting with weak turbulence.
    %         Physics of Fluids 33, 086111 (featured article). DOI: 10.1063/5.0059948.
    %
    %     [2] Cuadra, A., Vera, M., Di Renzo, M., & Huete, C. (2023). Linear Theory
    %         of Hypersonic Shocks Interacting with Turbulence in Air. In 2023 AIAA
    %         SciTech Forum, National Harbor, USA. DOI: 10.2514/6.2023-0075.
    %
    %     [3] https://in.mathworks.com/matlabcentral/fileexchange/28135-standard-atmosphere-functions?requestedDomain=

    % Default values
    T = 300;              % Temperature of the reactants [K]
    p = 1.01325;          % Pressure of the reactants [bar]
    species = {'N2', 'O2', 'Ar', 'CO2'}; % Reactants
    moles = [3.7276, 1.0000, 0.0447, 0.0015]; % Molar composition of the reactants
    LS = 'air_ions';      % List of possible products
    Npoints = 5e3;        % Approximate number of points for the calculations
    MAX_MACH = 30;        % Maximum pre-shock Mach number for the calculations
    FLAG_ATMOS = false;   % Flag indicating to compute pre-shock state for a given altitude [m]
    FLAG_RESULTS = false; % Flag to show results in the command window

    % Unpack inputs
    for i = 1:2:nargin
        switch lower(varargin{i})
            case {'t', 'tr'}
                T = varargin{i + 1};
            case {'p', 'pr'}
                p = varargin{i + 1};
            case {'species', 'reactants'}
                species = varargin{i + 1};
            case {'moles', 'n'}
                moles = varargin{i + 1};
            case {'ls', 'list_products', 'list products'}
                LS = varargin{i + 1};
            case {'z', 'altitude'}
                z = varargin{i + 1};
                FLAG_ATMOS = true;
            case {'flag_results'}
                FLAG_RESULTS = varargin{i + 1};
            case {'npoints', 'points'}
                Npoints = varargin{i + 1};
            case {'mach', 'max_mach'}
                MAX_MACH = varargin{i + 1};
        end

    end

    % Compute pre-shock conditions based on altitude above sea level [m]
    if FLAG_ATMOS
        [~, ~, T, p] = atmos(z(i));
        % Change pressure units
        p = convert_Pa_to_bar(p);
    end

    % Initialization
    self = App(LS);
    self.Misc.FLAG_RESULTS = FLAG_RESULTS;

    % Set initial conditions
    self = set_prop(self, 'TR', T, 'pR', p);
    self.PD.S_Oxidizer = species;
    self.PD.N_Oxidizer = moles;

    % Additional inputs
    initial_velocity_sound = compute_sound(T, p, self.PD.S_Oxidizer, self.PD.N_Oxidizer, 'self', self);
    u1 = [logspace(2, 4, 4 / 3 * Npoints), logspace(4, 5, round(Npoints / 100))];
    u1(u1 > MAX_MACH * initial_velocity_sound) = [];
    u1(u1 < initial_velocity_sound) = [];  
    self = set_prop(self, 'u1', u1);

    % Solve problem
    self = solve_problem(self, 'SHOCK_I');

    % Get jump conditions
    [R, P, T, M1, M2, Gammas] = get_parameters(self);
end

% SUB-PASS FUNCTIONS
function [R, P, T, M1, M2, Gammas] = get_parameters(self, varargin)
    % Get jump conditions

    % Default values
    FLAG_COMPOSITION = false;
    % Unpack 
    if nargin > 1
        FLAG_COMPOSITION = varargin{1};
    end

    % Compute jump conditions
    R = cell2vector(self.PS.strP, 'rho') ./ cell2vector(self.PS.strR, 'rho');
    P = cell2vector(self.PS.strP, 'p') ./ cell2vector(self.PS.strR, 'p');
    T = cell2vector(self.PS.strP, 'T') ./ cell2vector(self.PS.strR, 'T');
    M1 = cell2vector(self.PS.strR, 'u') ./ cell2vector(self.PS.strR, 'sound');
    M2 = cell2vector(self.PS.strP, 'v_shock') ./ cell2vector(self.PS.strP, 'sound');
    Gammas = compute_Gammas_frozen(M1, R, P);

    % Remove last value to have vectors of the same size
    R = R(1:end-1);
    P = P(1:end-1);
    T = T(1:end-1);
    M1 = M1(1:end-1);
    M2 = M2(1:end-1);

    if ~FLAG_COMPOSITION
        return
    end

    ind = find_ind(self.S.LS, {'O2','N2','O','N','NO','eminus','Nplus','Oplus'});
    Xi = cell2vector(self.PS.strP, 'Xi');

    Xi = Xi(:, 1:end-1);
    XO2 = Xi(ind(1), :);
    XN2 = Xi(ind(2), :);
    XO = Xi(ind(3), :);
    XN = Xi(ind(4), :);
    XNO = Xi(ind(5), :);
    Xeminus = Xi(ind(6), :);
    XNplus = Xi(ind(7), :);
    XOplus = Xi(ind(8), :);

    dXO2dM1 = compute_first_derivative(XO2, M1);
    dXN2dM1 = compute_first_derivative(XN2, M1);
    dXOdM1 = compute_first_derivative(XO, M1);
    dXNdM1 = compute_first_derivative(XN, M1);
    dXNOdM1 = compute_first_derivative(XNO, M1);
    dXeminusdM1 = compute_first_derivative(Xeminus, M1);
    dXNplusdM1 = compute_first_derivative(XNplus, M1);
    dXOplusdM1 = compute_first_derivative(XOplus, M1);

    M1_new = M1(1:end-1);
    
    xx = linspace(M1_new(1), M1_new(end), 250);
    dXO2dM1 = interp1(M1_new, dXO2dM1, xx);
    dXN2dM1 = interp1(M1_new, dXN2dM1, xx);
    dXOdM1 = interp1(M1_new, dXOdM1, xx);
    dXNdM1 = interp1(M1_new, dXNdM1, xx);
    dXNOdM1 = interp1(M1_new, dXNOdM1, xx);
    dXeminusdM1 = interp1(M1_new, dXeminusdM1, xx);
    dXNplusdM1 = interp1(M1_new, dXNplusdM1, xx);
    dXOplusdM1 = interp1(M1_new, dXOplusdM1, xx);

    M1dXO2dM1 = [M1_new; dXO2dM1]';
    M1dXN2dM1 = [M1_new; dXN2dM1]';
    M1dXOdM1 = [M1_new; dXOdM1]';
    M1dXNdM1 = [M1_new; dXNdM1]';
    M1dXNOdM1 = [M1_new; dXNOdM1]';
    M1dXeminusdM1 = [M1_new; dXeminusdM1]';
    M1dXNplusdM1 = [M1_new; dXNplusdM1]';
    M1dXOplusdM1 = [M1_new; dXOplusdM1]';

    M1_new = xx;

end