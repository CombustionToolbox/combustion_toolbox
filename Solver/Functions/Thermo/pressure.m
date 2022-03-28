function value = pressure(mix)
    % Get the pressure [bar] in the mixture
    %
    % Args:
    %     mix (struct):  Properties of the mixture
    %
    % Returns:
    %     value (float): Pressure [bar] in the mixture

    value = mix.p;
end