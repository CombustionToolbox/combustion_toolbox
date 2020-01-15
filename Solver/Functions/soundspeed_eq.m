function sound = soundspeed_eq(str0,phi,P0,T0,E,S,C,M,PD,TN,strThProp)
PD.ProblemType = 'TP';

T2 = 1.02*T0;
T1 = 0.98*T0;

str1 = SolveProblemTP_TV(str0,phi,P0,T1,E,S,C,M,PD,TN,strThProp);
s1 = str1.S*1e3;

str2 = SolveProblemTP_TV(str0,phi,P0,T2,E,S,C,M,PD,TN,strThProp);
s2 = str2.S*1e3;

DSDT = (s2 - s1)/(T2 - T1);

P1 = 0.98*P0;
P2 = 1.02*P0;

str1 = SolveProblemTP_TV(str0,phi,P1,T0,E,S,C,M,PD,TN,strThProp);
s1 = str1.S*1e3;

str2 = SolveProblemTP_TV(str0,phi,P2,T0,E,S,C,M,PD,TN,strThProp);
s2 = str2.S*1e3;

DSDP = (s2 - s1)/(P2 - P1)*1e-5;

DTDP = -DSDP/DSDT;

TA = T0 + DTDP*(P1-P0)*1e5;

strA = SolveProblemTP_TV(str0,phi,P1,TA,E,S,C,M,PD,TN,strThProp);
rhoA = strA.rho;

TB = T0 + DTDP*(P2-P0)*1e5;
strB = SolveProblemTP_TV(str0,phi,P2,TB,E,S,C,M,PD,TN,strThProp);
rhoB = strB.rho;

DRHODP = (P2-P1)*1e5/(rhoB-rhoA);

sound = sqrt(DRHODP);

