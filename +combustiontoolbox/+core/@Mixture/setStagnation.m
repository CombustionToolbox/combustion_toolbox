function mix = setStagnation(mix, varargin)
    % Update mixture properties to stagnation conditions
    %
    % Args:
    %     mix (Mixture):  Mixture object
    %
    % Returns:
    %     mix (Mixture):  Mixture object at stagnation conditions
    %
    % Example:
    %     mix = setStagnation(mix)

    % Import packages
    import combustiontoolbox.common.*
    import combustiontoolbox.equilibrium.*
    
    % Unpack inputs
    if nargin > 1
        FLAG_THERMO_CHEM = varargin{1};
    else
        FLAG_THERMO_CHEM = ~mix.FLAG_REACTION;
    end

    % Shallow copy
    mix = mix.copy;

    % Thermochemical frozen case
    if FLAG_THERMO_CHEM
        % Additional shallow copy
        mix1 = mix.copy;

        solver = EquilibriumSolver('problemType', 'TP', 'FLAG_TCHEM_FROZEN', true, 'FLAG_RESULTS', false);
        mix.T = mix.T * (1 + 0.5 * (mix.gamma - 1) * mix.mach^2);
        mix.p = mix.p * (1 + 0.5 * (mix.gamma - 1) * mix.mach^2)^(mix.gamma / (mix.gamma - 1));
        solver.equilibrateT(mix1, mix, mix.T);

        % Set flow velocity [m/s]
        mix.u = 0; mix.uShock = 0; mix.uNormal = 0;
        return
    end

    % Definitions
    s_target = mix.s;
    h_target = mix.h + 0.5 * mix.u^2 * mix.mi;
    solver = EquilibriumSolver('problemType', 'HP', 'FLAG_TCHEM_FROZEN', false, 'FLAG_RESULTS', false);

    % Set initial guess for stagnation pressure (assuming calorically perfect gas) [bar] 
    p = mix.p * (1 + 0.5 * (mix.gamma_s - 1) * mix.mach^2)^(mix.gamma_s / (mix.gamma_s - 1));

    % Set flow velocity [m/s]
    mix.u = 0; mix.uShock = 0; mix.uNormal = 0;

    % Get equilibrium state assuming an isentropic process
    STOP = 1; it = 0;
    delta = 1;
    mix0 = mix;

    while STOP > solver.tol0 && it < solver.itMax
        % Update iteration
        it = it + 1;
        
        % Define state
        mix = mix0.copy();
        mix.p = p;
        mix.h = h_target;

        % Solve equilibrium state for the current pressure
        solver.solve(mix);
        [s, p] = getState(mix);

        % Solve equilibrium state for the perturbed pressure
        p_perturb = mix.p * 1.01;
        mix.h = h_target;
        mix.p = p_perturb;
        solver.solve(mix);
        [s_perturb, p_perturb] = getState(mix);
    
        % Compute residual and first derivative
        f0 = s - s_target;
        df = (s_perturb - s) / (p_perturb - p);
        
        % Update pressure
        p = abs(p - delta * f0 / df);
        
        % Compute stop criteria
        STOP = max(abs(f0 / s_target), abs(f0));
    end

    solver.solve(mix);

    % Check if maximum iterations were reached
    if it >= solver.itMax
        warning('Maximum iterations reached without full convergence.');
    end

end

% SUB-PASS FUNCTIONS
function [entropy, pressure] = getState(mix)
    % Get state
    entropy = mix.s;
    pressure = mix.p;
end