function [x, STOP] = steff(self, mix1, pP, field, x0)
    % Find the temperature [K] (root) for the set chemical transformation at equilibrium using the Steffenson-Aitken method
    %
    % Args:
    %     self (struct):  Data of the mixture, conditions, and databases
    %     mix1 (struct):  Properties of the initial mixture
    %     pP (float):     Pressure [bar]
    %     field (str):    Fieldname in Problem Description (PD)
    %     x0 (float):     Guess temperature [K]
    % Returns:
    %     - x (float):    Temperature at equilibrium [K]
    %     - STOP (float): Relative error [-] 

    if any(strcmpi(self.PD.ProblemType, {'TP', 'TV'}))
        x = get_transformation(self, 'TP');
        STOP = 0;
        return
    end

    it = 0; g = 1.0; STOP = 1.0;
    
    while STOP > self.TN.tol0 && it < self.TN.itMax
        it = it + 1;
        [g, g_rel]= get_gpoint(self, mix1, pP, field, x0);
        fx = abs(g - x0);
        g_aux  = get_gpoint(self, mix1, pP, field, fx);
        fx2 = abs(g_aux - fx);
        if abs(fx2 - 2*fx + x0) > self.TN.tol0
            x = get_point_aitken(x0, [fx, fx2]);
        else
            x = fx;
        end

        STOP = max(abs((x - fx) / x), abs(g_rel));
        x0 = x;
    end
    print_error_root(it, self.TN.itMax, x, STOP);
end