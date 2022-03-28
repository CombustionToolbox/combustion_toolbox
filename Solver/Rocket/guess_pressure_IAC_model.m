function pressure = guess_pressure_IAC_model(mix)
    % Compute pressure guess [bar] at the throat considering an Infinite-Area-Chamber (IAC) 
    %
    % This method is based on Gordon, S., & McBride, B. J. (1994). NASA reference publication,
    % 1311.
    %
    % Args:
    %     mix (struct): Properties of the mixture
    %
    % Returns:
    %     pressure (float): Pressure [bar]

    pressure = mix.p * ((mix.gamma_s + 1) / 2) ^ -(mix.gamma_s / (mix.gamma_s - 1));
end