function sound = soundspeed_eq(app, str0, phi, P0, T0)
app.PD.ProblemType = 'TP';

T2 = 1.02*T0;
T1 = 0.98*T0;

str1 = SolveProblemTP_TV(app, str0, phi, P0, T1);
s1 = str1.S*1e3;

str2 = SolveProblemTP_TV(app, str0, phi, P0, T2);
s2 = str2.S*1e3;

DSDT = (s2 - s1)/(T2 - T1);

P1 = 0.98*P0;
P2 = 1.02*P0;

str1 = SolveProblemTP_TV(app, str0, phi, P1, T0);
s1 = str1.S*1e3;

str2 = SolveProblemTP_TV(app, str0, phi, P2, T0);
s2 = str2.S*1e3;

DSDP = (s2 - s1)/(P2 - P1)*1e-5;

DTDP = -DSDP/DSDT;

TA = T0 + DTDP*(P1-P0)*1e5;

strA = SolveProblemTP_TV(app, str0, phi, P1, TA);
rhoA = strA.rho;

TB = T0 + DTDP*(P2-P0)*1e5;
strB = SolveProblemTP_TV(app, str0, phi, P2, TB);
rhoB = strB.rho;

DRHODP = (P2-P1)*1e5/(rhoB-rhoA);

sound = sqrt(DRHODP);

