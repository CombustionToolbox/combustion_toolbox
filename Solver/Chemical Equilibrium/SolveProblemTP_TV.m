function [strP] = SolveProblemTP_TV(strR,phi,pP,TP,E,S,C,M,PD,TN,strThProp)
% CALCULATE EQUILIBRIUM AT DEFINED T AND P (TP)
%                       OR
% CALCULATE EQUILIBRIUM AT DEFINED T AND CONSTANT V (TV)
% INPUT:
%   strR  = Prop. of reactives (phi,species,...)
%   phi   = velocity upstream       [-]
%   pP    = pressure of products    [bar]
%   TP    = temperature of products [K]
% OUTPUT:
%   strP  = Prop. of products (phi,species,...)

[NCO2P_0,NCOP_0,NH2OP_0,NH2P_0,NO2P_0,NN2P_0,NHeP_0,NArP_0,NCgrP_0,phi_c,FLAG_SOOT] = CalculateProductsCC(strR.NatomE,phi,TP,E.Elements,TN.factor_c,PD.Fuel,strThProp);
P = SetSpecies(C.M0.Value,S.List_Compute_Species,[NCO2P_0,NCOP_0,NH2OP_0,NH2P_0,NO2P_0,NN2P_0,NHeP_0,NArP_0,NCgrP_0]',TP,S.idx_fixed,strThProp);
if strcmpi(PD.CompleteOrIncomplete,'INCOMPLETE')
    % Compute number of moles of M.minor_products
    P = CalculateProductsIC(P,TP,pP,strR.v,phi,M.minor_products,phi_c,FLAG_SOOT,C.M0.Value,C.A0.Value,E,S,C,M,PD,TN,strThProp);
    % Compute properties of all species
    P = SetSpecies(C.M0.Value,S.List_Compute_Species,P(S.idx_all,1),TP,S.idx_all,strThProp);
end

if strfind(PD.ProblemType,'P') == 2 % PD.ProblemType = 'TP', 'HP', or 'SP'
    strP = ComputeProperties(C.A0.Value,P,pP,TP,E.ind_C,E.ind_H);
elseif strfind(PD.ProblemType,'V') == 2 % PD.ProblemType = 'TV', 'EV', or 'SV'
    NP = sum(P(:,1).*(1-P(:,10)));
    pP = (NP*TP*8.3144598/(strR.v/1000))/1e5;
    strP = ComputeProperties(C.A0.Value,P,pP,TP,E.ind_C,E.ind_H);
end
end