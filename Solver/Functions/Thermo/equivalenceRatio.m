function value = equivalenceRatio(mix)
    % Get the equivalence ratio of the initial mixture [-]
    value = mix.phi;
    if ischar(value)
        value = nan;
    end
end