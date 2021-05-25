function [N_IC, DeltaNP] = Equilibrium(app, N_CC, phi, pP, TP, vR)
% Generalized Gibbs minimization method                             
N0 = app.C.N0.Value;
A0 = app.C.A0.Value;
R0TP = app.C.R0 * TP; % [J/(mol)]
% Initialization
NatomE = dot(N_CC, A0);
NP_0 = 0.1;
NP = NP_0;
it = 0;
% itMax = 500
itMax = 50 + round(S.NS/2);
SIZE = -log(app.C.tolN);
e = 0.;
STOP = 1.;
end