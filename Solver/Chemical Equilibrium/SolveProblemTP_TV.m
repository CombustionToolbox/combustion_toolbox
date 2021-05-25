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
phi_c0 = 2/(PD.Fuel.x-PD.Fuel.z)*(PD.Fuel.x+PD.Fuel.y/4-PD.Fuel.z/2);
[NCO2P_0,NCOP_0,NH2OP_0,NH2P_0,NO2P_0,NN2P_0,NHeP_0,NArP_0,NCgrP_0,phi_c0,FLAG_SOOT] = CalculateProductsCC(strR.NatomE,phi,phi_c0,TP,pP,E.elements,TN.factor_c,PD.Fuel,strThProp);
P = SetSpecies(C.M0.value,S.LS,[NCO2P_0,NCOP_0,NH2OP_0,NH2P_0,NO2P_0,NN2P_0,NHeP_0,NArP_0,NCgrP_0]',TP,S.ind_fixed,strThProp);
if strcmpi(PD.CompleteOrIncomplete,'INCOMPLETE')
    % Compute number of moles of M.minor_products
                                             
    [P,DeltaNP] = CalculateProductsIC(P,TP,pP,strR.v,phi,M.minors_products,phi_c0,FLAG_SOOT,C.M0.value,C.A0.value,E,S,C,M,PD,TN,strThProp);
    % Compute properties of all species
    P = SetSpecies(C.M0.value,S.LS,P(S.ind_all,1),TP,S.ind_all,strThProp);
else
    DeltaNP = 0;
end

if strfind(PD.ProblemType,'P') == 2 % PD.ProblemType = 'TP', 'HP', or 'SP'
    strP = ComputeProperties(C.A0.value,P,pP,TP,E.ind_C,E.ind_H);
elseif strfind(PD.ProblemType,'V') == 2 % PD.ProblemType = 'TV', 'EV', or 'SV'
    NP = sum(P(:,1).*(1-P(:,10)));
    pP = (NP*TP*8.3144598/(strR.v/1000))/1e5;
    strP = ComputeProperties(C.A0.value,P,pP,TP,E.ind_C,E.ind_H);
end
strP.phi_c = phi_c0;
strP.error_moles = DeltaNP;
end