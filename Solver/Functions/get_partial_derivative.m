function value = get_partial_derivative(self, mix)
    % Get value of the partial derivative for the set problem type [kJ/K] (HP, EV) or [kJ/K^2] (SP, SV)
    %
    % Args:
    %     self (struct): Data of the mixture, conditions, and databases
    %     mix (struct):  Properties of the mixture
    %
    % Returns:
    %     value (float): Value of the partial derivative for the set problem type [kJ/K] (HP, EV) or [kJ/K^2] (SP, SV)

    if strcmpi(self.PD.ProblemType, 'HP')
        value = mix.cP;
    elseif strcmpi(self.PD.ProblemType, 'EV')
        value = mix.cV;
    elseif strcmpi(self.PD.ProblemType, 'SP')
        value = mix.cP / mix.T;
    elseif strcmpi(self.PD.ProblemType, 'SV')
        value = mix.cV / mix.T;
    end
    value = value * 1e-3; % [kJ/K] (HP, EV) or [kJ/K^2] (SP, SV)
end