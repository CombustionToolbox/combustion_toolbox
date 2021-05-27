function [strP] = SolveProblemSV(app, strR, phi, vP)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE ISENTROPIC COMPOSITION AT CONSTANT V
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
%   strR  = Prop. of reactives (phi,species,...)
%   phi   = equivalence ratio      [-]
%   vP    = volumen of products    [m^3]
% OUTPUT:
%   strP  = Prop. of products (phi,species,...)

% constants;
% if C.firstrow, TP = 0.75*strR.v/vP*strR.T; end
% DeltaT = 1;
% 
% strR.v = vP;
% pP = strR.p;
% 
% while abs(DeltaT) > 1e-2
%     
% %     relaxT = exp(-((TP-1000)/2500).^2);
%     relaxT = 1;
%         
%     strP = SolveProblemTP_TV(strR,phi,pP,TP,E,S,C,M,PD,TN,strThProp);
% 
%     DS  = strP.S - strR.S;
%     cVP = strP.cV;
% 
%     DeltaT = - DS*TP/(cVP);
% 
%     TP = TP + relaxT*DeltaT;
%     
%     if TP > 6000, TP = 6000; end
% 
% end

if vP>strR.v % expansion
    TP_l = 250;
    TP_r = strR.T;
else % compression
    TP_l = strR.T;
    TP_r = strR.T+1000;
end
strR.v = vP;
pP = strR.p;

if app.C.firstrow
    
    strP = SolveProblemTP_TV(app, strR, phi, pP ,TP_l);
    if isnan(strP.S)
        TP_l = TP_l+100;
        strP = SolveProblemTP_TV(app, strR, phi, pP, TP_l);
    end
    f_l  = strP.S - strR.S;
    strP = SolveProblemTP_TV(app, strR, phi, pP, TP_r);
    f_r  = strP.S - strR.S;
    
    if f_l*f_r > 0 || (isnan(f_l) && isnan(f_r))
        TP = strR.T+1500;
    elseif abs(f_l)<abs(f_r) || abs(f_l)>=abs(f_r)
        TP = TP_r - (TP_r-TP_l)/(f_r-f_l)*f_r;
        strP = SolveProblemTP_TV(app, strR, phi, pP, TP);
        f  = strP.S - strR.S;
        %           fun = griddedInterpolant([f_l f f_r],[TP_l TP TP_r],'makima');
        %           TP = interp1([f_l f_r f],[TP_l TP_r TP],0);
        TP = interp1([f_l f_r],[TP_l TP_r],0);
        %           TP = fun(0);
    elseif isnan(f_l) && ~isnan(f_r)
        TP = TP_r-100;
    elseif ~isnan(f_l) && isnan(f_r)
        TP = TP_l+100;
    else
        TP = TP_r - (TP_r-TP_l)/(f_r-f_l)*f_r;
        strP = SolveProblemTP_TV(app, strR, phi, pP, TP);
        f  = strP.S - strR.S;
        %           fun = griddedInterpolant([f_l f f_r],[TP_l TP TP_r],'makima');
        TP = interp1([f_l f_r f],[TP_l TP_r TP],0);
        %           TP = fun(0);
    end
end
DeltaT = 1;
tol0 = 1e-10;
it = 0;
% f = 1;
itMax = 100;
while (abs(DeltaT) > 1e-3 || abs(f) > 1e-3) && it<itMax
    it = it+1;
    strP = SolveProblemTP_TV(app, strR, phi, pP, TP);
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