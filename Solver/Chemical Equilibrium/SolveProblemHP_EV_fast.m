function [strP] = SolveProblemHP_EV_fast(strR,phi,pP,strP,E,S,C,M,PD,TN,strThProp)
% CALCULATE ADIABATIC T AND COMPOSITION AT CONSTANT P (HP)
%                           OR
% CALCULATE EQUILIBRIUM COMPOSITION AT ADIABATIC T AND CONSTANT V (EV)
% INPUT:
%   strR  = Prop. of reactives (phi,species,...)
%   phi   = equivalence ratio       [-]
%   pP    = pressure of products    [bar]
% OUTPUT:
%   strP  = Prop. of products (phi,species,...)
DeltaT = 1;
tol0 = 1e-10;
it = 0;
Q = 1;
itMax = 30;
TP = strP.T;
if strcmp(PD.ProblemType,'HP')
    while (abs(DeltaT) > 1e-2 || abs(Q) > 1e-2) && it<itMax
        it = it+1;
        strP = SolveProblemTP_TV(strR,phi,pP,TP,E,S,C,M,PD,TN,strThProp);
        Q  = strP.h - strR.h;
        gx = abs(Q-TP);
        strP_aux = SolveProblemTP_TV(strR,phi,pP,gx,E,S,C,M,PD,TN,strThProp);
        Q_aux  = strP_aux.h - strR.h;
        gx2 = abs(Q_aux-gx);
        if abs(gx2-2*gx+TP) > tol0
            TP = TP - (gx-TP)^2/(gx2-2*gx+TP);
        else
            TP = gx;
        end
        DeltaT = abs(Q_aux-Q)/(1 + abs(Q_aux));
    end
else
    while (abs(DeltaT) > 1e-2 || abs(Q) > 1e-2) && it<itMax
        it = it+1;
        strP = SolveProblemTP_TV(strR,phi,pP,TP,E,S,C,M,PD,TN,strThProp);
        Q  = strP.e - strR.e;
        gx = abs(Q-TP);
        strP_aux = SolveProblemTP_TV(strR,phi,pP,gx,E,S,C,M,PD,TN,strThProp);
        Q_aux  = strP_aux.e - strR.e;
        gx2 = abs(Q_aux-gx);
        if abs(gx2-2*gx+TP) > tol0
            TP = TP - (gx-TP)^2/(gx2-2*gx+TP);
        else
            TP = gx;
        end
        DeltaT = abs(Q_aux-Q)/(1 + abs(Q_aux));
    end
end

strP.error_problem = max(abs(DeltaT),abs(Q));

if it>=itMax
    fprintf('****************************\n');
    fprintf('** Solution not converged **\n');
    fprintf('** phi   =  %4.2f         **\n',phi);
    fprintf('** Temp  =  %4.2f         **\n',TP);
    fprintf('** Error =  %4.2f%%         **\n',abs(DeltaT)*100);
    fprintf('** It    =  %4.d          **\n',it);
    fprintf('****************************\n');
end
% fprintf('** It    =  %4.d          **\n',it);