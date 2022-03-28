function [gpoint, gpoint_relative] = get_gpoint(self, mix1, pP, field, x0)
    % Get fixed point of a function based on the chemical transformation
    %
    % Args:
    %     self (struct):  Data of the mixture, conditions, and databases
    %     mix1 (struct):  Properties of the initial mixture
    %     pP (float):     Pressure [bar]
    %     field (str):    Fieldname in Problem Description (PD)
    %     x0 (float):     Guess temperature [K]
    % Returns:
    %     gpoint (float): Fixed point of the function [kJ] (HP, EV) or [kJ/K] (SP, SV)
    %     gpoint_relative (float): Fixed relative point of the function [kJ] (HP, EV) or [kJ/K] (SP, SV)

    mix2 = equilibrate_T(self, mix1, pP, x0);
    gpoint = (mix2.(field) - mix1.(field));
    gpoint_relative = gpoint / (mix2.(field));
    if strcmpi(field, 's')
        gpoint = gpoint * 1e3;
    end
end