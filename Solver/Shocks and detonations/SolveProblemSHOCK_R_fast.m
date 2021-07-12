function [str1, str2, str3] = SolveProblemSHOCK_R_fast(app, strR, phi, p1, T1, u1, str2, str3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE PLANAR INCIDENT SHOCK STATE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[str1, str2] = SolveProblemSHOCK_I_fast(app, strR, phi, p1, T1, u1, str2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE POST-REFLECTED SHOCK STATE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
%   R  = Prop. Species upstram (phi,species,...)
%   p1 = pressure state 1    [bar]
%   T1 = temperature state 1 [K]
%   u1 = velocity state 1    [m/s]
% OUTPUT:
%   str1 = Prop. Species state 1
%   str2 = Prop. Species state 2
%   str3 = Prop. Species state 3
% NOTES:
%   Index 1: state 1, upstream
%   Index 2: state 2, incident
%   Index 3: state 3, reflected
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% help SolveProblemSHOCK_I
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GLOBAL CONFIGURATION OF shock_incident.m
% TOLERANCES
TN.ERRFT = 1e-4;
TN.ERRFV = 1e-4;
% INITIAL ASSUMPTIONS OF ERROR
deltaT = 1000; 
deltaV = 1000; 
% MISCELANEOUS
j = 0;     % initilize looping variable
nfrec = 0; % frequency of sampling results
pP = p1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD VARIABLES
V1 = str1.V;
T1 = str1.T;
u1 = str1.u;
r1 = str1.rho;
h1 = str1.h/str1.mi*1e3;

V2 = str2.V;
T2 = str2.T;
u2 = str2.u;
p2 = str2.p;
r2 = str2.rho;
h2 = str2.h/str2.mi*1e3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
strR = str2; % propertiers mixture at state 2
UI = u2;     % velocity incident shock [m/s]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PRELIMINARY GUESS
u2 = sqrt((p2*1e5-p1*1e5)*(V1-V2)); % particle velocity [m/s]
w2 = UI - u2; % velocity in shock fixed frame [m/s]

V = 1/str3.rho;
r = 1/V;
T = str3.T;
h = str3.h/str3.mi*1e3;
p = str3.p;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START LOOP
while((abs(deltaT) > TN.ERRFT*T) || (abs(deltaV) > TN.ERRFV*V))
    j = j + 1;
    if(j == 500)
        disp(['shock_reflected did not converge for u1 = ',num2str(u1)])
        return
    end       
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %CALCULATE FH & FP FOR GUESS 1
    FH = h -h2- 0.5*u2^2* (r/r2+1)/(r/r2-1);
    FP = p*1e5 - p2*1e5 - r*u2^2/(r/r2-1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %TEMPERATURE PERTURBATION
    DT = T*0.02;
    Tper = T + DT;
    Vper = V;
    rper = 1/Vper;
    
%     stateper;
    strP = state(app, strR, rper, Tper, phi, pP);
    hper = strP.h/strP.mi*1e3;
    pper = strP.p;
    
    %CALCULATE FHX & FPX FOR "IO" STATE
    FHX = hper -h2- 0.5*u2^2* (rper/r2+1)/(rper/r2-1);
    FPX = pper*1e5 - p2*1e5 - rper*u2^2/(rper/r2-1);

    %ELEMENTS OF JACOBIAN
    DFHDT = (FHX-FH)/DT;
    DFPDT = (FPX-FP)/DT;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %VOLUME PERTURBATION
    DV = 0.02*V;
    Vper = V + DV;
    Tper = T;
    rper = 1/Vper;
    
%     stateper;
    strP = state(app, strR, rper, Tper, phi, pP);
    hper = strP.h/strP.mi*1e3;
    pper = strP.p;
    
    %CALCULATE FHX & FPX FOR "IO" STATE
    FHX = hper -h2- 0.5*u2^2*(rper/r2+1)/(rper/r2-1);
    FPX = pper*1e5 - p2*1e5 - rper*u2^2/(rper/r2-1);
   
    %ELEMENTS OF JACOBIAN
    DFHDV = (FHX-FH)/DV;
    DFPDV = (FPX-FP)/DV;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %USE MATLAB MATRIX INVERTER
    J = [DFHDT DFHDV; DFPDT DFPDV];
    a = [-FH; -FP];
    b = J\a;
    deltaT = b(1);
    deltaV = b(2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %CHECK & LIMIT CHANGE VALUES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %TEMPERATURE
    DTM = 0.2*T;
    if (abs(deltaT) > DTM)
        deltaT = DTM*deltaT/abs(deltaT);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %VOLUME
    V3X = V + deltaV;
    if (V3X > V2)
        DVM = 0.5*(V2 - V);
    else
        DVM = 0.2*V;
    end
    if (abs(deltaV) > DVM)
        deltaV = DVM*deltaV/abs(deltaV);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %MAKE THE CHANGES
    T = T + deltaT;
    V = V + deltaV;
    r = 1/V;
%     state;
    strP = state(app, strR, r, T, phi, pP);
    h = strP.h/strP.mi*1e3;
    p = strP.p;
%     UR = (p-p2)*1e5/(u2*r2)-u2; % NASA
    UR = (p-p2)*1e5/(u2*r2)-u2; % CALTECH same as Gaseq
%     UR = sqrt(((p-p2)*1e5)/((r/r2-1)*r)); % Notes, same as Gaseq
    % PRINT SAMPLING RESULTS
%     if(mod(j,nfrec)==0)
%         disp('----------------------------------------------------')
%         disp('Solution Planar Reflected Shock')
%         fprintf('  T = %8.3f [K] \n',T)
%         fprintf('  p = %8.3f [bar] \n',p)
%         fprintf('  u = %8.3f [m/s] \n',UR)
%         fprintf('  r = %8.3f [kg/m3] \n',r)
%         fprintf('  h = %8.3f [kJ/kg] \n',h*1e-3)      
%     end
end
% SOLUTION PLANAR INCIDENT SHOCK 
T3 = T; % temperature downstream     [K]
V3 = V; % specific volume downstream [m3/kg]
p3 = p; % pressure downstream        [bar]
h3 = h; % enthalpy downstream        [J/kg]
r3 = r; % density downstream         [kg/m3]
u3 = UR; % velocity downstream       [m/s]
% PRINT RESULTS
% disp('-----------------------------------------------------------')
% fprintf('Number of iterations = %d \n',j)
% fprintf('maxerror = %.2e\n',max(abs(deltaV),abs(deltaT)))
% disp('-----------------------------------------------------------')
% fprintf('\t\t\t STATE 1\t\tSTATE 2\t\t\tSTATE 3\n')
% fprintf('\t\t\t INITIAL\t\tINCIDENT\t\tREFLECTED\n')
% fprintf('T [K]    \t %6.3f\t\t%6.3f \t\t%6.3f\n',T1,T2,T3)
% fprintf('p [bar]\t %11.3f\t\t%7.3f \t\t%7.3f\n',p1,p2,p3)
% fprintf('u [m/s]\t %11.3f\t\t%7.3f \t\t%7.3f\n',u1,UI,UR)
% fprintf('r [kg/m3]\t %7.3f\t\t%7.3f \t\t%7.3f\n',r1,r2,r3)
% fprintf('h [kJ/kg]\t %7.3f\t\t%7.3f \t\t%7.3f\n',h1*1e-3,h2*1e-3,h3*1e-3)
% disp('-----------------------------------------------------------')
% SAVE STATE
str3 = strP;
str3.u = u3;
str3.T = T3;
str3.V = V3;
str3.p = p3;
str3.h = h3*str3.mi/1e3;

str3.error_problem = max(abs(deltaT),abs(deltaV));
end

function strP = state(app, strR, r, T, phi, pP)
% Calculate frozen state given T & rho
strR.v = strR.mi/r*1e3;
TP = T; % vP = vR (computed from R);
% Equilibrium composition at defined T and constant v
app.PD.ProblemType = 'TV';
strP = SolveProblemTP_TV(app, strR, phi, pP, TP);
end