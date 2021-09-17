function strP = state(self, strR, r, T, pP)
% Calculate frozen state given T & rho
strR.v = strR.mi/r*1e3;
self.PD.TP.value = T;
% Equilibrium composition at defined T and constant v
self.PD.ProblemType = 'TV';
strP = equilibrate(self, strR, pP);
end