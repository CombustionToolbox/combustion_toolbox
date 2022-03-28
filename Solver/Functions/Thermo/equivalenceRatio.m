function value = equivalenceRatio(mix)
    % Get the equivalence ratio of the initial mixture [-]
    %
    % Args:
    %     mix (struct):  Properties of the mixture
    %
    % Returns:
    %     value (float): Equivalence ratio of the initial mixture [-]

    value = mix.phi;
    if ischar(value)
        value = nan;
    end
end