function value = soundspeed(mix)
    % Get the speed of sound [m/s] in the mixture
    %
    % Args:
    %     mix (struct):  Properties of the mixture
    %
    % Returns:
    %     value (float): Speed of sound [m/s] in the mixture

    value = mix.sound;
    if isnan(val)
        value = Inf;
    end
end