function value = velocity_relative(mix)
    % Get the velocity of the gases relative to the shock front [m/s] in the mixture
    %
    % Args:
    %     mix (struct):  Properties of the mixture
    %
    % Returns:
    %     value (float): Velocity of the gases relative to the shock front [m/s] in the mixture

    value = mix.u;
end