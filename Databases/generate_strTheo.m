function strTheo = generate_strTheo()
% Database for theoretical computation of the jump conditions of a diatomic
% species only considering dissociation.

% N2
strTheo.N2.Tr = 2.87;
strTheo.N2.Tv = 3390;
strTheo.N2.Td = 113000;
strTheo.N2.G  = 4^2;
strTheo.N2.m  = 2.3259*1e-26;
% O2
strTheo.O2.Tr = 2.08;
strTheo.O2.Tv = 2270;
strTheo.O2.Td = 59500;
strTheo.O2.G  = 5^2/3;
strTheo.O2.m  = 2.6567*10^-26;
% H2
strTheo.H2.Tr = 87.53;
strTheo.H2.Tv = 6338;
strTheo.H2.Td = 51973;
strTheo.H2.G  = 2^2;
strTheo.H2.m  = 0.16735*10^-26;
% I2
strTheo.I2.Tr = 0.0538;
strTheo.I2.Tv = 308;
strTheo.I2.Td = 17897;
strTheo.I2.G  = 4^2;
strTheo.I2.m  = 21.072*10^-26;
% F2
strTheo.F2.Tr = 1.27;
strTheo.F2.Tv = 1320;
strTheo.F2.Td = 18633;
strTheo.F2.G  = 4^2;
strTheo.F2.m  = 3.1548*10^-26;
% Cl2
strTheo.Cl2.Tr = 0.0346;
strTheo.Cl2.Tv = 2270;
strTheo.Cl2.Td = 28770;
strTheo.Cl2.G  = 4^2;
strTheo.Cl2.m  = 5.8871*10^-26;