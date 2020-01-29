function [strP] = SolveProblemEV(strR,phi,pP,E,S,C,M,PD,TN,strThProp)
% CALCULATE EQUILIBRIUM COMPOSITION AT ADIABATIC T AND CONSTANT V (EV)
%
% INPUT:
%   strR  = Prop. of reactives (phi,species,...)
%   phi   = equivalence ratio       [-]
%   pP    = pressure of products    [bar]
% OUTPUT:
%   strP  = Prop. of products (phi,species,...)

if C.firstrow
    TP_l = 800;
    TP_r = 1500;
    strP = SolveProblemTP_TV(strR,phi,pP,TP_l,E,S,C,M,PD,TN,strThProp);
    if isnan(strP.e)
        TP_l = TP_l+100;
        strP = SolveProblemTP_TV(strR,phi,pP,TP_l,E,S,C,M,PD,TN,strThProp);
    end
    Q_l  = strP.e - strR.e;
    strP = SolveProblemTP_TV(strR,phi,pP,TP_r,E,S,C,M,PD,TN,strThProp);
    Q_r  = strP.e - strR.e;
    
    if Q_l*Q_r > 0 || (isnan(Q_l) && isnan(Q_r))
        TP = 2500;
    elseif abs(Q_l)<abs(Q_r) || abs(Q_l)>=abs(Q_r)
        TP = TP_r - (TP_r-TP_l)/(Q_r-Q_l)*Q_r;
        strP = SolveProblemTP_TV(strR,phi,pP,TP,E,S,C,M,PD,TN,strThProp);
        Q  = strP.e - strR.e;
        %           fun = griddedInterpolant([Q_l Q Q_r],[TP_l TP TP_r],'makima');
        %           TP = interp1([Q_l Q_r Q],[TP_l TP_r TP],0);
        TP = interp1([Q_l Q_r],[TP_l TP_r],0);
        %           TP = fun(0);
    elseif isnan(Q_l) && ~isnan(Q_r)
        TP = TP_r-100;
    elseif ~isnan(Q_l) && isnan(Q_r)
        TP = TP_l+100;
    else
        TP = TP_r - (TP_r-TP_l)/(Q_r-Q_l)*Q_r;
        strP = SolveProblemTP_TV(strR,phi,pP,TP,E,S,C,M,PD,TN,strThProp);
        Q  = strP.e - strR.e;
        %           fun = griddedInterpolant([Q_l Q Q_r],[TP_l TP TP_r],'makima');
        TP = interp1([Q_l Q_r Q],[TP_l TP_r TP],0);
        %           TP = fun(0);
    end
end

DeltaT = 1;
tol0 = 1e-10;
% Q = 1;
itMax = 30;
it = 0;
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