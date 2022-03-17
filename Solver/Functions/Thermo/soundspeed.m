function val = soundspeed(mix)
    % Get the speed of sound [m/s]
    val = mix.sound;
    if isnan(val)
        val = Inf;
    end
end