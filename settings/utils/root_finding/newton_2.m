function [x, STOP] = newton_2(f, fprime, x0)
    % Find the temperature [K] (root) for the set chemical transformation at equilibrium using the Newton-Raphson method
    %
    % Args:
    %     self (struct): Data of the mixture, conditions, and databases
    %     mix1 (struct): Properties of the initial mixture
    %     pP (float): Pressure [bar]
    %     field (str): Fieldname in Problem Description (PD)
    %     x0 (float): Guess temperature [K]
    %     guess_moles (float): Guess moles final mixture
    %
    % Returns:
    %     Tuple containing
    %
    %     - x (float): Temperature at equilibrium [K]
    %     - STOP (float): Relative error [-] 
    %     - guess_moles (struct): Guess moles final mixture

    it = 0; STOP = 1.0;

    while STOP > 1e-4 && it < 30
        it = it + 1;
        
        x = abs(x0 - f0(x0) / fprime0(x0));

        STOP = abs((x - x0) / x);
        x0 = x;
    end
end