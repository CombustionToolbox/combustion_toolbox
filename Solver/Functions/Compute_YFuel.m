function Yi_Fuel = Compute_YFuel(mix, mix_Fuel)
    % Compute fuel mass fraction [-]
    %
    % Args:
    %     mix (struct):       Properties of the mixture (fuel + oxidizer + inerts)
    %     mix_Fuel (struct):  Properties of the mixture (fuel)
    %
    % Returns:
    %     Yi_Fuel (float):    Mass fractions of the fuel mixture

    for i=length(mix):-1:1
        Yi_Fuel(i) = mix_Fuel.mi/ mix{i}.mi;
    end
end
