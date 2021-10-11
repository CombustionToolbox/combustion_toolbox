function DB_Theo = generate_DB_Theo()
% Database for theoretical computation of the jump conditions of a diatomic
% species only considering dissociation.

% N2
DB_Theo.N2.Tr = 2.87;
DB_Theo.N2.Tv = 3390;
DB_Theo.N2.Td = 113000;
DB_Theo.N2.G  = 4^2;
DB_Theo.N2.m  = 2.3259*1e-26;
% O2
DB_Theo.O2.Tr = 2.08;
DB_Theo.O2.Tv = 2270;
DB_Theo.O2.Td = 59500;
DB_Theo.O2.G  = 5^2/3;
DB_Theo.O2.m  = 2.6567*10^-26;
% H2
DB_Theo.H2.Tr = 87.53;
DB_Theo.H2.Tv = 6338;
DB_Theo.H2.Td = 51973;
DB_Theo.H2.G  = 2^2;
DB_Theo.H2.m  = 0.16735*10^-26;
% I2
DB_Theo.I2.Tr = 0.0538;
DB_Theo.I2.Tv = 308;
DB_Theo.I2.Td = 17897;
DB_Theo.I2.G  = 4^2;
DB_Theo.I2.m  = 21.072*10^-26;
% F2
DB_Theo.F2.Tr = 1.27;
DB_Theo.F2.Tv = 1320;
DB_Theo.F2.Td = 18633;
DB_Theo.F2.G  = 4^2;
DB_Theo.F2.m  = 3.1548*10^-26;
% Cl2
DB_Theo.Cl2.Tr = 0.0346;
DB_Theo.Cl2.Tv = 2270;
DB_Theo.Cl2.Td = 28770;
DB_Theo.Cl2.G  = 4^2;
DB_Theo.Cl2.m  = 5.8871*10^-26;