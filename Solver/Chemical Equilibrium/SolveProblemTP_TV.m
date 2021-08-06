function [strP] = SolveProblemTP_TV(self, strR, pP, TP)
% CALCULATE EQUILIBRIUM AT DEFINED T AND P (TP)
%                       OR
% CALCULATE EQUILIBRIUM AT DEFINED T AND CONSTANT V (TV)
% INPUT:
%   strR  = Prop. of reactives (phi,species,...)
%   phi   = Equivalence ratio [-]
%   TP    = Products temperature [K]
% OUTPUT:
%   strP  = Prop. of products (phi,species,...)

% Compute number of moles 
[N, STOP] = Equilibrium(self, pP, TP, strR);
% Compute properties of all species
P = SetSpecies(self, self.S.LS, N(:, 1), TP);

if strfind(self.PD.ProblemType, 'P') == 2 % PD.ProblemType = 'TP', 'HP', or 'SP'
    strP = ComputeProperties(self, P, pP, TP);
else % PD.ProblemType = 'TV', 'EV', or 'SV'
    NP = sum(P(:,1).*(1-P(:,10)));
    pP = (NP*TP*8.3144598/(strR.v/1000))/1e5;
    strP = ComputeProperties(self, P, pP, TP);
end
strP.phi_c = Compute_phi_c(self.PD.Fuel);
strP.error_moles = STOP;
end
