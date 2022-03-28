function value = temperature(mix)
    % Get the temperature [K] in the mixture
    %
    % Args:
    %     mix (struct):  Properties of the mixture
    %
    % Returns:
    %     value (float): Temperature [K] in the mixture

    value = mix.T;
end