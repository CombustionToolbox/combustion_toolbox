function pressure = guess_pressure_IAC_model(mix)
    % Compute pressure guess [bar] for Infinite-Area-Chamber (IAC) 
    pressure = mix.p * ((mix.gamma_s + 1) / 2) ^ -(mix.gamma_s / (mix.gamma - 1));
end