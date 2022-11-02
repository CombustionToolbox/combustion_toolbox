function sound = soundspeed_eq(self, mix, P0, T0)
    % Compute speed of sound at equilibrium
    %
    % Args:
    %     self (struct): Data of the mixture, conditions, and databases
    %     mix (struct): Struct mix with all the properties of the mixture
    %     P0 (float): Pressure [bar]
    %     T0 (float): Temperature [K]
    %
    % Returns:
    %     sound (float): sound speed [m/s]

    self.PD.ProblemType = 'TP';

    T2 = 1.02 * T0;
    T1 = 0.98 * T0;

    str1 = equilibrate_T(self, mix, P0, T1, []);
    guess_moles = str1.Xi * str1.N;
    s1 = str1.S * 1e3;

    str2 = equilibrate_T(self, mix, P0, T2, guess_moles);
    s2 = str2.S * 1e3;

    DSDT = (s2 - s1) / (T2 - T1);

    P1 = 0.98 * P0;
    P2 = 1.02 * P0;

    str1 = equilibrate_T(self, mix, P1, T0, guess_moles);
    s1 = str1.S * 1e3;

    str2 = equilibrate_T(self, mix, P2, T0, guess_moles);
    s2 = str2.S * 1e3;

    DSDP = (s2 - s1) / (P2 - P1) * 1e-5;

    DTDP = -DSDP / DSDT;

    TA = T0 + DTDP * (P1 - P0) * 1e5;

    strA = equilibrate_T(self, mix, P1, TA, guess_moles);
    rhoA = strA.rho;

    TB = T0 + DTDP * (P2 - P0) * 1e5;
    strB = equilibrate_T(self, mix, P2, TB, guess_moles);
    rhoB = strB.rho;

    DRHODP = (P2 - P1) * 1e5 / (rhoB - rhoA);

    sound = sqrt(DRHODP);
end
