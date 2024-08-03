function pressure = rocketGuessThroatIAC(mix)
    % Compute pressure guess [bar] at the throat considering an Infinite-Area-Chamber (IAC)
    %
    % Args:
    %     mix (Mixture): Properties of the mixture
    %
    % Returns:
    %     pressure (float): Pressure at the throat [bar]
    %
    % Example:
    %     pressure = rocketGuessThroatIAC(mix)

    pressure = mix.p / ((mix.gamma_s + 1) / 2)^(mix.gamma_s / (mix.gamma_s - 1));
end