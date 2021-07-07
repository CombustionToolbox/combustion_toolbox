function [strP] = SolveProblemTP_TV(app, strR, phi, pP, TP)
% CALCULATE EQUILIBRIUM AT DEFINED T AND P (TP)
%                       OR
% CALCULATE EQUILIBRIUM AT DEFINED T AND CONSTANT V (TV)
% INPUT:
%   strR  = Prop. of reactives (phi,species,...)
%   phi   = Equivalence ratio [-]
%   pP    = Products pressure [bar]
%   TP    = Products temperature [K]
% OUTPUT:
%   strP  = Prop. of products (phi,species,...)

% Abbreviations ---------------------
E = app.E;
S = app.S;
C = app.C;
PD = app.PD;
% -----------------------------------

if strcmpi(app.PD.solver, 'SEGREGATED')
    [N_CC, phi_c0, FLAG_SOOT, STOP] = CalculateProductsCC(app, strR, phi, pP, TP);
    P = SetSpecies(C.M0.value, S.LS, N_CC', TP, S.ind_fixed, app.strThProp);
    if strcmpi(PD.CompleteOrIncomplete,'INCOMPLETE')
    % N_CC matrix with number of moles and swtCondesated of each species
    N_CC = P(:, [1,10]);
    % Compute number of moles of M.minor_products
    [N_IC, STOP] = CalculateProductsIC(app, N_CC, phi, pP, TP, strR.v, phi_c0, FLAG_SOOT);
    end
else
    if strcmpi(app.PD.solver, 'GIBBS')
        [N_IC, STOP] = Equilibrium(app, pP, TP, strR);
%         [N_IC, TP, STOP] = Equilibrium_HP(app, pP, strR); % TEST
    elseif strcmpi(app.PD.solver, 'GIBBS SOOT')
        [N_IC, STOP] = Equilibrium_soot(app, pP, TP, strR);
    elseif strcmpi(app.PD.solver, 'GIBBS REDUCED')
    [N_IC, STOP] = Equilibrium_reduced(app, pP, TP, strR);
    end
    phi_c0 = Compute_phi_c(app.PD.Fuel);
end

if strcmpi(PD.CompleteOrIncomplete,'COMPLETE') && strcmpi(app.PD.solver, 'SEGREGATED')
    STOP = 0;
else
    % Compute properties of all species
    P = SetSpecies(C.M0.value, S.LS, N_IC(S.ind_all, 1), TP, S.ind_all, app.strThProp);
end

if strfind(PD.ProblemType,'P') == 2 % PD.ProblemType = 'TP', 'HP', or 'SP'
    strP = ComputeProperties(C.A0.value,P,pP,TP,E.ind_C,E.ind_H);
elseif strfind(PD.ProblemType,'V') == 2 % PD.ProblemType = 'TV', 'EV', or 'SV'
    NP = sum(P(:,1).*(1-P(:,10)));
    pP = (NP*TP*8.3144598/(strR.v/1000))/1e5;
    strP = ComputeProperties(C.A0.value,P,pP,TP,E.ind_C,E.ind_H);
end
strP.phi_c = phi_c0;
strP.error_moles = STOP;
end
