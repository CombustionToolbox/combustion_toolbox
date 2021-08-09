function strP = Equilibrate(varargin)
    try
        self = varargin{1};
        strR = varargin{2}{1,1};
        pP = varargin{3};
        if nargin == 4, strP = varargin{4}; else, strP = []; end
        % get attribute of the specified transformations
        attr_name = get_attr_name(self);
        % compute initial guess
        guess = get_guess(self, strR, pP, attr_name, strP);
        % root finding: find the value x that satisfies f(x) = strP.xx(x) - strR.xx = 0
        [x, ERR] = root_finding(self, strR, pP, attr_name, guess);
        % compute properties
        strP = equilibrate_T(self, strR, pP, x);
        strP.error_problem = ERR;
    catch
        error("An exception occurred: error Equilibrate.m")
    end
end


%%% NESTED FUNCTIONS
function strP = equilibrate_T(self, strR, pP, TP)
    % Compute number of moles 
    [N, DeltaNP] = Equilibrium(self, pP, TP, strR);
    % Compute properties of all species
    P = SetSpecies(self, self.S.LS, N(:, 1), TP);

    if strfind(self.PD.ProblemType, 'P') == 2
        strP = ComputeProperties(self, P, pP, TP);
    else
        NP = sum(P(:, 1) * (1 - P(:, 10)));
        pP = (NP * TP * self.C.R0 / (strR.v/1e3)) / 1e5;
        strP = ComputeProperties(self, P, pP, TP);
    end    
    strP.error_moles = DeltaNP;
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


function guess = get_guess(self, strR, pP, attr_name, strP)
    if any(strcmpi(self.PD.ProblemType, {'TP', 'TV'}))
        guess = get_transformation(self, 'TP');
    elseif ~isempty(strP)
        guess = strP.T;
    else
        guess = steff_guess(self, strR, pP, attr_name);
    end
end


function [x, ERR] = root_finding(self, strR, pP, attr_name, x0)
    [x, ERR] = self.TN.root_method(self, strR, pP, attr_name, x0);
end
        