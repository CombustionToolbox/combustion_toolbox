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

    mix1 = get_struct(mix1);
    [mix2, guess_moles] = unpack(varargin);
    % get attribute xx of the specified transformations
    attr_name = get_attr_name(self);
    % compute initial guess
    guess = get_guess(self, mix1, pP, attr_name, mix2);
    % root finding: find the value x that satisfies f(x) = mix2.xx(x) - mix1.xx = 0
    [T, ERR] = root_finding(self, mix1, pP, attr_name, guess);
    % compute properties
    mix2 = equilibrate_T(self, mix1, pP, T, guess_moles);
    % check convergence in case the problemType is TP (defined Temperature and Pressure)
    print_convergence(mix2.error_moles, self.TN.tolN, mix2.error_moles_ions, self.TN.tol_pi_e, self.PD.ProblemType)
    % save error - root finding
    mix2.error_problem = ERR;
end

%%% SUB-PASS FUNCTIONS
function [mix2, guess_moles] = unpack(value)
    if ~isempty(value)
        mix2 = get_struct(value{1});
        guess_moles = mix2.Xi * mix2.N;
    else
        mix2 = [];
        guess_moles = [];
    end
end

function str = get_struct(var)
    try
        str = var{1,1};
    catch
        str = var;
    end
end

function attr_name = get_attr_name(self)
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

function guess = get_guess(self, mix1, pP, attr_name, mix2)
    if any(strcmpi(self.PD.ProblemType, {'TP', 'TV'}))
        guess = get_transformation(self, 'TP');
    elseif ~isempty(mix2)
        guess = mix2.T;
    else
        guess = steff_guess(self, mix1, pP, attr_name);
    end
end

function [x, ERR] = root_finding(self, mix1, pP, attr_name, x0)
    [x, ERR] = self.TN.root_method(self, mix1, pP, attr_name, x0);
end

function print_convergence(STOP, TOL, STOP_ions, TOL_ions, ProblemType)
    if strcmpi(ProblemType, 'TP')
        if STOP > TOL
            fprintf('***********************************************************\n')
            fprintf('Convergence error number of moles:   %.2f\n', STOP);
        end
        if STOP_ions > TOL_ions
            fprintf('***********************************************************\n')
            fprintf('Convergence error in charge balance: %.2f\n', STOP_ions);
        end
    end
end
        