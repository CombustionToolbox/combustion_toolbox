function Z = compressibility_factor(mix)
    % Compute compressibility factor of the mixture [-]
    %
    % Args:
    %     mix (struct): Properties of the mixture
    %
    % Returns:
    %     Z (float): Compressibility factor [-]

    Z = pressure(mix) * volume(mix) / (molesGas(mix) * 8.31446261815324 * temperature(mix)) * 1e5;
end
