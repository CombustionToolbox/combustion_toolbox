function W = Compute_W(mix, mix_Fuel)
    % Compute non-dimensional slope of the variation of the density over the equivalence ratio [-]
    % Paper: Cuadra (2020). Effect of equivalence ratio fluctuations on planar detonation discontinuities. J. Fluid Mech 903, A30.
    %
    % Args:
    %     mix (struct):       Properties of the mixture (fuel + oxidizer + inerts)
    %     mix_Fuel (struct):  Properties of the mixture (fuel)
    %
    % Returns:
    %     W (float):          Non-dimensional slope of the variation of the density over the equivalence ratio [-]

    for i=length(mix):-1:1
        Y_Fuel = mix_Fuel.mi / mix{i}.mi;
        W_Air  = mix{i}.mi - mix_Fuel.mi;
        W_Fuel = mix_Fuel.mi;
        W(i)  = (1 - (W_Air / W_Fuel)) / (1 - (1 - W_Air / W_Fuel) * Y_Fuel);    
    end
end
