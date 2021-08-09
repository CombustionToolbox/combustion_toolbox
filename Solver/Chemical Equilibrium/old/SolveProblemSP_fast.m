function [strP] = SolveProblemSP_fast(app, strR, phi, pP, strP)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE ISENTROPIC COMPOSITION AT CONSTANT P
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
%   strR  = Prop. of reactives (phi,species,...)
%   phi   = equivalence ratio       [-]
%   pP    = pressure of products    [bar]
% OUTPUT:
%   strP  = Prop. of products (phi,species,...)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% help SolveProblemSP_test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% constants;
% if C.firstrow, TP = 2000; end
% 
% DeltaT = 1;
% 
% while abs(DeltaT) > 1e-2
%     
%     relaxT = 0.85*exp(-((TP-1000)/2500).^2);
%         
%     strP = SolveProblemTP_TV(strR,phi,pP,TP,E,S,C,M,PD,TN,strThProp);
% 
%     DS  = strP.S - strR.S;
%     cPP = strP.cP;
% 
%     DeltaT = - DS*TP/(cPP);
% 
%     TP = TP + relaxT*DeltaT;
%     
%     if TP > 6000, TP = 6000; end
% 
% end
DeltaT = 1;
tol0 = 1e-10;
it = 0;
f = 1;
itMax = 20;
TP = strP.T;
while (abs(DeltaT) > 1e-3 || abs(f) > 1e-3) && it<itMax
    it = it+1;
    strP = SolveProblemTP_TV(app, strR, phi, pP, TP);
    if isnan(strP.S)
        TP = 1.05*TP;
        strP = SolveProblemTP_TV(app, strR, phi, pP, TP);
        if isnan(strP.S)
           TP = 0.9*TP;
            strP = SolveProblemTP_TV(app, strR, phi, pP, TP); 
        end
    end
    f  = strP.S - strR.S;
    gx = abs(f-TP);
    strP_aux = SolveProblemTP_TV(app, strR, phi, pP, gx);
    f_aux  = strP_aux.S - strR.S;
    gx2 = abs(f_aux-gx);
    if abs(gx2-2*gx+TP) > tol0
        TP = TP - (gx-TP)^2/(gx2-2*gx+TP);
    else
        TP = gx;
    end
    DeltaT = abs(f_aux-f)/(1 + abs(f_aux));
end

strP.error_problem = max(abs(DeltaT),abs(f));

if it>=itMax
    fprintf('****************************\n');
	fprintf('** Solution not converged **\n');
    fprintf('** phi   =  %4.2f         **\n',phi);
    fprintf('** Temp  =  %4.2f         **\n',TP);
    fprintf('** Error =  %4.2f%%         **\n',abs(DeltaT)*100);
    fprintf('** It    =  %4.d          **\n',it);
    fprintf('****************************\n');
end