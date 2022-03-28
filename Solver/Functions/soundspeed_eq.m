function sound = soundspeed_eq(self, mix, phi, P0, T0)
    % Compute speed of sound at equilibrium
    %
    % Args:
    %     self (struct):  Data of the mixture, conditions, and databases
    %     mix (struct):   Struct mix with all the properties of the mixture
    %     phi (float):    Equivalence ratio [-]
    %     P0 (float):     Pressure [bar]
    %     T0 (float):     Temperature [K]
    %
    % Returns:
    %     sound (float):  sound speed [m/s]
    
    self.PD.ProblemType = 'TP';

    T2 = 1.02*T0;
    T1 = 0.98*T0;

    str1 = SolveProblemTP_TV(self, mix, phi, P0, T1);
    s1 = str1.S*1e3;

    str2 = SolveProblemTP_TV(self, mix, phi, P0, T2);
    s2 = str2.S*1e3;

    DSDT = (s2 - s1)/(T2 - T1);

    P1 = 0.98*P0;
    P2 = 1.02*P0;

    str1 = SolveProblemTP_TV(self, mix, phi, P1, T0);
    s1 = str1.S*1e3;

    str2 = SolveProblemTP_TV(self, mix, phi, P2, T0);
    s2 = str2.S*1e3;

    DSDP = (s2 - s1)/(P2 - P1)*1e-5;

    DTDP = -DSDP/DSDT;

    TA = T0 + DTDP*(P1-P0)*1e5;

    strA = SolveProblemTP_TV(self, mix, phi, P1, TA);
    rhoA = strA.rho;

    TB = T0 + DTDP*(P2-P0)*1e5;
    strB = SolveProblemTP_TV(self, mix, phi, P2, TB);
    rhoB = strB.rho;

    DRHODP = (P2-P1)*1e5/(rhoB-rhoA);

    sound = sqrt(DRHODP);
end

