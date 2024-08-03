function mix2 = equilibrate(obj, mix2, varargin)
    % Obtain properties at equilibrium for the given thermochemical transformation
    %
    % Args:
    %     obj (EquilibriumSolver): EquilibriumSolver object
    %     mix2 (Mixture): Mixture considering a thermochemical frozen gas
    %
    % Optional Args:
    %     * mix2Guess (Mixture): Mixture object of a previous calculation
    %
    % Returns:
    %     mix2 (Mixture): Mixture at chemical equilibrium for the given thermochemical transformation
    %
    % Examples:
    %     * mix2 = equilibrate(EquilibriumSolver(), mix2)
    %     * mix2 = equilibrate(EquilibriumSolver(), mix2, mix2Guess)
    
    % Initialization
    mix1 = mix2.copy();

    if obj.FLAG_TCHEM_FROZEN
        % TO BE IMPLEMENTED
        mix2.errorProblem = 0;
        return
    end
    
    % Default
    mixGuess = [];

    % Unpack TGuess
    if nargin > 2
        mixGuess = varargin{1};
    end

    % Get attribute xx of the specified transformations
    attributeName = getAttribute(obj);

    % Get temperature and chemical composition TGuess
    [TGuess, molesGuess] = getGuess(obj, mix1, mix2, mixGuess, attributeName);

    % Root finding: find the value x that satisfies f(x) = mix2.xx(x) - mix1.xx = 0
    [T, STOP, molesGuess] = rootFinding(obj, mix1, mix2, attributeName, TGuess, molesGuess);

    % Compute properties
    equilibrateT(obj, mix1, mix2, T, molesGuess);

    % Check convergence in case the problemType is TP (defined Temperature and Pressure)
    print_convergence(mix2.errorMoles, obj.tolGibbs, mix2.errorMolesIons, obj.tolMultiplierIons, obj.problemType)

    % Save error from root finding algorithm
    mix2.errorProblem = STOP;
end

% SUB-PASS FUNCTIONS
function attributeName = getAttribute(obj)
    % Get attribute of the problem type
    switch upper(obj.problemType)
        case {'TP', 'TV'}
            attributeName = 'T';
        case 'HP'
            attributeName = 'h';
        case 'EV'
            attributeName = 'e';
        case {'SP', 'SV'}
            attributeName = 's';
    end

end

function [TGuess, molesGuess] = getGuess(obj, mix1, mix2, mixGuess, attributeName)
    % Get initial estimates for temperature and molar composition

    % Initialization
    if any(strcmpi(obj.problemType, {'TP', 'TV'}))
        TGuess = mix2.T;
        molesGuess = [];

        if ~isempty(mixGuess)
            molesGuess = mixGuess.Xi * mixGuess.N;
        end
        
        return
    end

    if ~isempty(mixGuess)
        TGuess = mixGuess.T;
        molesGuess = mixGuess.Xi * mixGuess.N;
    else
        TGuess = obj.regulaGuess(mix1, mix2, attributeName);
        molesGuess = [];
    end

end

function [x, STOP, molesGuess] = rootFinding(obj, mix1, mix2, attributeName, x0, molesGuess)
    % Calculate the temperature value that satisfied the problem conditions
    % using the @rootMethod
    [x, STOP, molesGuess] = obj.rootMethod(obj, mix1, mix2, attributeName, x0, molesGuess);
end

function print_convergence(STOP, TOL, STOP_ions, TOL_ions, ProblemType)
    % Print tolerance error if the convergence criteria was not satisfied

    if ~strcmpi(ProblemType, 'TP')
        return
    end

    if STOP > TOL
        fprintf('***********************************************************\n')
        fprintf('Convergence error number of moles:   %1.2e\n', STOP);
    end

    if STOP_ions > TOL_ions
        fprintf('***********************************************************\n')
        fprintf('Convergence error in charge balance: %1.2e\n', STOP_ions);
    end

end

function mix = set_volume_SV(obj, mix)
    % If the problem type is SV, the product's volume is based on the given v_P/v_R ratio
    if ~strcmpi(obj.problemType, 'SV')
        return
    end

    %mix.v = mix.v * obj.PD.vP_vR.value;
    fprintf('\nto be clarified\n')
end