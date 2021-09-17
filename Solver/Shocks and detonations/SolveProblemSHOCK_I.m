function [str1,str2] = SolveProblemSHOCK_I(app, strR, p1, T1, u1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE PLANAR INCIDENT SHOCK WAVE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
%   strR  = Prop. Species upstram (phi,species,...)
%   p1    = pressure upstream    [bar]
%   T1    = temperature upstream [K]
%   u1    = velocity upstream    [m/s]
% OUTPUT:
%   str1  = Prop. Species upstream
%   str2  = Prop. Species downstream
% NOTES:
%   Index 1: state 1, upstream
%   Index 2: state 2, downstream
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% help SolveProblemSHOCK_I
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Abbreviations ---------------------
TN = app.TN;
% -----------------------------------
% INITIAL ASSUMPTIONS OF ERROR
deltaT = 1000; 
deltaV = 1000;
% MISCELANEOUS
j = 0;     % initilize looping variable
nfrec = 0; % frequency of sampling results
pP = p1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE PROPERTIES OF THE REACTIVE AT THE GIVEN CONDITIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assuming perfect gases
r1 = strR.rho;           % density upstream  [kg/m3]
V1 = 1/r1;               % specific volume   [m3/kg]
h1 = strR.h/strR.mi*1e3; % enthalpy upstream [J/kg]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PRELIMINARY GUESS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VOLUME BOUND RATIO
if isfield(TN,'guess_shock')
    TN.volumeBoundRation = V1/TN.guess_shock(3);
else
    TN.volumeBoundRation = 5; 
end
V = V1/TN.volumeBoundRation;   % specifiv volume downstream [m3/kg]
r = 1/V;                       % density downstream         [kg/m3] 
p = p1*1e5 + r1*u1^2*(1-V/V1); % pressure downstream        [Pa]
p = p/1e5;                     % pressure downstream        [bar]
T = T1*p*V/(p1*V1);            % Temperature downstream     [K]
% COMPUTE PROPERTIES OF THE PRODUCTS AT THE GIVEN CONDITIONS
% state;
strP = state(app, strR, r, T, pP);
h = strP.h/strP.mi*1e3; % enthalpy upstream   [J/kg]
p = strP.p;             % pressure downstream [bar]
u = u1*r1/r;            % velocity downstream [m/s]. Continuity equation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while((abs(deltaT) > TN.ERRFT*T) || (abs(deltaV) > TN.ERRFV*V)) 
    j = j + 1;
    if(j == 500)
        disp(['shock_incident did not converge for u1 = ',num2str(u1)])
        return
    end           
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CALCULATE FH & FP FOR GUESS 1
    FH = (h+0.5*u^2)-(h1+0.5*u1^2);
    FP = (p*1e5+r*u^2)-(p1*1e5+r1*u1^2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %TEMPERATURE PERTURBATION
    DT = 0.02*T;
    Tper = T+DT;
    Vper = V;
    rper = 1/Vper;

    strP = state(app, strR, rper, Tper, pP);
    uper = u1*r1/rper;
    hper = strP.h/strP.mi*1e3;
    pper = strP.p;
    % CALCULATE FHX & FPX FOR "IO" STATE
    FHX = (hper+0.5*uper^2)-(h1+0.5*u1^2);
    FPX = (pper*1e5+rper*uper^2)-(p1*1e5+r1*u1^2);
    %ELEMENTS OF JACOBIAN
    DFHDT = (FHX-FH)/DT;
    DFPDT = (FPX-FP)/DT;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % VOLUME PERTURBATION
    DV = 0.02*V;
    Tper = T;
    Vper = V + DV;
    rper = 1/Vper;

    strP = state(app, strR, rper, Tper, pP);
    uper = u1*r1/rper;
    hper = strP.h/strP.mi*1e3;
    pper = strP.p;
    % CALCULATE FHX & FPX FOR "IO" STATE
    FHX = (hper+0.5*uper^2)-(h1+0.5*u1^2);
    FPX = (pper*1e5+rper*uper^2)-(p1*1e5+r1*u1^2);
    % ELEMENTS OF JACOBIAN
    DFHDV = (FHX-FH)/DV;
    DFPDV = (FPX-FP)/DV;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % USE MATLAB MATRIX INVERTER
    J = [DFHDT DFHDV; DFPDT DFPDV];
    a = [-FH; -FP];
    b = J\a;
    deltaT = b(1);
    deltaV = b(2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CHECK & LIMIT CHANGE VALUES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TEMPERATURE
    DTM = 0.2*T;
    if (abs(deltaT) > DTM)
        deltaT = DTM*deltaT/abs(deltaT);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % VOLUME
    V2X = V + deltaV;
    if (V2X > V1)
        DVM = 0.5*(V1 - V);
    else
        DVM = 0.2*V;
    end
    if (abs(deltaV) > DVM)
        deltaV = DVM*deltaV/abs(deltaV);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MAKE THE CHANGES
    T = T + deltaT;
    V = V + deltaV;
    r = 1/V;
    
    strP = state(app, strR, r, T, pP);
    h = strP.h/strP.mi*1e3;
    p = strP.p;
    u = u1*r1/r;
    % PRINT SAMPLING RESULTS
%     if(mod(j,nfrec)==0)
%         disp('----------------------------------------------------')
%         disp('Solution Planar Incident Shock')
%         fprintf('it %d, maxerror = %.2e:\n',j,max(abs(deltaV),abs(deltaT)))
%         fprintf('  T = %8.3f [K] \n',T)
%         fprintf('  p = %8.3f [bar] \n',p)
%         fprintf('  u = %8.3f [m/s] \n',u)
%         fprintf('  r = %8.3f [kg/m3] \n',r)
%         fprintf('  h = %8.3f [kJ/kg] \n',h*1e-3)        
%     end
end
% SOLUTION PLANAR INCIDENT SHOCK 
T2 = T; % temperature downstream     [K]
V2 = V; % specific volume downstream [m3/kg]
p2 = p; % pressure downstream        [bar]
h2 = h; % enthalpy downstream        [J/kg]
r2 = r; % density downstream         [kg/m3]
u2 = u; % velocity downstream        [m/s]
%%%%%
% u2 = sqrt((p2*1e5-p1*1e5)*(V1-V2))
%%%
% PRINT RESULTS
% disp('-----------------------------------------------------------')
% fprintf('Number of iterations = %d \n',j)
% fprintf('maxerror = %.2e\n',max(abs(deltaV),abs(deltaT)))
% disp('-----------------------------------------------------------')
% fprintf('\t\t\t STATE 1\t\tSTATE 2\n')
% fprintf('\t\t\t INITIAL\t\tINCIDENT\n')
% fprintf('T [K]    \t %6.3f\t\t%6.3f\n',T1,T2)
% fprintf('p [bar]\t %11.3f\t\t%7.3f\n',p1,p2)
% fprintf('u [m/s]\t %11.3f\t\t%7.3f\n',u1,u2)
% fprintf('r [kg/m3]\t %7.3f\t\t%7.3f\n',r1,r2)
% fprintf('h [kJ/kg]\t %7.3f\t\t%7.3f\n',h1*1e-3,h2*1e-3)
% disp('-----------------------------------------------------------')
% SAVE STATES
str1 = strR;
str1.u = u1;
str1.T = T1;
str1.V = V1;
str2 = strP;
str2.u = u2;
str2.T = T2;
str2.V = V2;
str2.p = p2;
str2.h = h2*str2.mi/1e3;

str2.error_problem = max(abs(deltaT),abs(deltaV));
end