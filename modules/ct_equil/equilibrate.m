function mix2 = equilibrate(self, mix1, pP, varargin)
    % Obtain properties at equilibrium for the set thermochemical transformation
    %
    % Args:
    %     self (struct): Data of the mixture, conditions, and databases
    %     mix1 (struct): Properties of the initial mixture
    %     pP (float):    Pressure [bar]
    %
    % Optional Args:
    %     mix2 (struct): Properties of the final mixture (previous calculation)
    %
    % Returns:
    %     mix2 (struct): Properties of the final mixture
    %
    % Example:
    %     mix2 = equilibrate(self, self.PS.strR{1}, 1.01325)

    mix1 = get_struct(mix1);
    mix2 = unpack(varargin);

    if self.TN.FLAG_TCHEM_FROZEN
        mix2 = mix1;
        mix2.p = pP;
        mix2.error_problem = 0;
        mix2.phi = get_phi(mix1);
        return
    end

    % get attribute xx of the specified transformations
    attr_name = get_attr_name(self);
    % compute initial guess
    [guess, guess_moles] = get_guess(self, mix1, pP, attr_name, mix2);
    % root finding: find the value x that satisfies f(x) = mix2.xx(x) - mix1.xx = 0
    [T, STOP, guess_moles] = root_finding(self, mix1, pP, attr_name, guess, guess_moles);
    % compute properties
    mix2 = equilibrate_T(self, mix1, pP, T, guess_moles);
    % check convergence in case the problemType is TP (defined Temperature and Pressure)
    print_convergence(mix2.error_moles, self.TN.tol_gibbs, mix2.error_moles_ions, self.TN.tol_pi_e, self.PD.ProblemType)
    % save error from root finding algorithm
    mix2.error_problem = STOP;
    % save equivalence ratio
    mix2.phi = get_phi(mix1);
end

%%% SUB-PASS FUNCTIONS
function mix2 = unpack(value)
    % Unpack input parameters
    if ~isempty(value)
        mix2 = get_struct(value{1});
    else
        mix2 = [];
    end

end

function str = get_struct(var)
    % Get struct
    try
        str = var{1, 1};
    catch
        str = var;
    end

end

function attr_name = get_attr_name(self)
    % Get attribute of the problem type
    if any(strcmpi(self.PD.ProblemType, {'TP', 'TV'}))
        attr_name = 'T';
    elseif any(strcmpi(self.PD.ProblemType, 'HP'))
        attr_name = 'h';
    elseif any(strcmpi(self.PD.ProblemType, 'EV'))
        attr_name = 'e';
    elseif any(strcmpi(self.PD.ProblemType, {'SP', 'SV'}))
        attr_name = 'S';
    end

end

function [guess, guess_moles] = get_guess(self, mix1, pP, attr_name, mix2)
    % Get initial estimates for temperature and molar composition
    if any(strcmpi(self.PD.ProblemType, {'TP', 'TV'}))
        guess = get_transformation(self, 'TP');
        guess_moles = [];
    elseif ~isempty(mix2)
        guess = mix2.T;
        guess_moles = mix2.Xi * mix2.N;
    else
        guess = regula_guess(self, mix1, pP, attr_name);
        guess_moles = [];
    end

end

function [x, STOP, guess_moles] = root_finding(self, mix1, pP, attr_name, x0, guess_moles)
    % Calculate the temperature value that satisfied the problem conditions
    % using the @root_method
    [x, STOP, guess_moles] = self.TN.root_method(self, mix1, pP, attr_name, x0, guess_moles);
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

function value = get_phi(mix1)
    % Get equivalence ratio
    try
        value = mix1.phi;
    catch
        value = [];
    end

end
