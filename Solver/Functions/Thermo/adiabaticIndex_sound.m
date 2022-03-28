function value = adiabaticIndex_sound(mix)
    % Get the adiabatic index [-] of the mixture from definition of sound velocity
    %
    % Args:
    %     mix (struct):  Properties of the mixture
    %
    % Returns:
    %     value (float): Adiabatic index [-] of the mixture

    value = mix.gamma_s;
end