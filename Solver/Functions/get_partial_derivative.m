function value = get_partial_derivative(str)
    if strcmpi(ProblemType, 'HP')
        value = str.cP;
    elseif strcmpi(ProblemType, 'EV')
        value = str.cV;
    elseif strcmpi(ProblemType, 'SP')
        value = str.cP / str.T;
    elseif strcmpi(ProblemType, 'SV')
        value = str.cV / str.T;
    end
    value = value * 1e-3; % units
end