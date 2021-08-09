function value = get_partial_derivative(self, str)
    if strcmpi(self.PD.ProblemType, 'HP')
        value = str.cP;
    elseif strcmpi(self.PD.ProblemType, 'EV')
        value = str.cV;
    elseif strcmpi(self.PD.ProblemType, 'SP')
        value = str.cP / str.T;
    elseif strcmpi(self.PD.ProblemType, 'SV')
        value = str.cV / str.T;
    end
    value = value * 1e-3; % units
end