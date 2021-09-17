function [str1,str2,guess] = SolveProblemDET_Jouget_inloop(strR,phi,p1,T1,volumeBoundRation,guess,E,S,C,M,PD,TN,strThProp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE CHAPMAN-JOUGET STATE (CJ UPPER STATE - STRONG DETONATION)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% help SolveProblemDET
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GLOBAL CONFIGURATION
% constants;
% TOLERANCES
TN.ERRFT = 1e-3;
TN.ERRFV = 1e-3;
% INITIAL ASSUMPTIONS OF ERROR
deltaT = 1000;
deltaV = 1000;
% VOLUME BOUND RATIO
% volumeBoundRation = 1.472;
% MISCELANEOUS
j = 0;     % initilize looping variable
nfrec = 0; % frequency of sampling results
pP = p1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assuming perfect gases
r1 = strR.rho;
V1 = 1/r1;
h1 = strR.h/strR.mi*1e3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PRELIMINARY GUESS
T = guess(1);
V = V1/volumeBoundRation;
r = 1/V;
% STATE
strP = state(strR,r,T,phi,pP,E,S,C,M,PD,TN,strThProp);
h = strP.h/strP.mi*1e3;
p = strP.p;
u = soundspeed_eq(strP,phi,p,T,E,S,C,M,PD,TN,strThProp);
u1=u*r/r1;
% u1 = guess(2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START LOOP
tic
while((abs(deltaT) > TN.ERRFT*T) || (abs(deltaV) > TN.ERRFV*V))
    j = j + 1;
    if(j == 500)
        disp('CJ_detonation did not converge')
        return
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %CALCULATE FH & FP FOR GUESS 1
    FH = (h+0.5*u^2)-(h1+0.5*u1^2);
    FP = (p*1e5+r*u^2)-(p1*1e5+r1*u1^2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %TEMPERATURE PERTURBATION
    DT = T*0.02;
    Tper = T + DT;
    Vper = V;
    rper = 1/Vper;   
    %     stateper;
    strP = state(strR,rper,Tper,phi,pP,E,S,C,M,PD,TN,strThProp);
    hper = strP.h/strP.mi*1e3;
    pper = strP.p;
    uper = soundspeed_eq(strP,phi,pper,Tper,E,S,C,M,PD,TN,strThProp);
%     u1=uper*r/r1;
    %CALCULATE FHX & FPX FOR "IO" STATE
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
    
    % 	stateper;
    strP = state(strR,rper,Tper,phi,pP,E,S,C,M,PD,TN,strThProp);
    hper = strP.h/strP.mi*1e3;
    pper = strP.p;
    uper = soundspeed_eq(strP,phi,pper,Tper,E,S,C,M,PD,TN,strThProp);
    %CALCULATE FHX & FPX FOR "IO" STATE
    FHX = (hper+0.5*uper^2)-(h1+0.5*u1^2);
    FPX = (pper*1e5+rper*uper^2)-(p1*1e5+r1*u1^2);
    %ELEMENTS OF JACOBIAN
    DFHDV = (FHX-FH)/DV;
    DFPDV = (FPX-FP)/DV;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %USE MATLAB MATRIX INVERTER
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
    % 	state;
    strP = state(strR,r,T,phi,pP,E,S,C,M,PD,TN,strThProp);
    h = strP.h/strP.mi*1e3;
    p = strP.p;
    u = soundspeed_eq(strP,phi,p,T,E,S,C,M,PD,TN,strThProp);
%     u1=u*r/r1;
    % PRINT SAMPLING RESULTS
    %     if(mod(j,nfrec)==0)
    %         disp('----------------------------------------------------')
    %         disp('Solution CJ Detonation')
    %         fprintf('it %d, maxerror = %.2e:\n',j,max(abs(deltaU),abs(deltaT)))
    %         fprintf('  T = %8.3f [K] \n',T)
    %         fprintf('  p = %8.3f [bar] \n',p)
    %         fprintf('  u = %8.3f [m/s] \n',u)
    %         fprintf('  u1 = %8.3f [m/s] \n',u1)
    %         fprintf('  r = %8.3f [kg/m3] \n',r)
    %         fprintf('  h = %8.3f [kJ/kg] \n',h*1e-3)
    %     end
end
j 
toc
% SOLUTION CJ DETONATION
T2 = T; % temperature downstream     [K]
V2 = V; % specific volume downstream [m3/kg]
p2 = p; % pressure downstream        [bar]
h2 = h; % enthalpy downstream        [J/kg]
r2 = r; % density downstream         [kg/m3]
u2 = u; % velocity downstream        [m/s]
% PRINT RESULTS
% disp('-----------------------------------------------------------')
% fprintf('Number of iterations = %d \n',j)
% fprintf('maxerror = %.2e\n',max(abs(deltaU),abs(deltaT)))
% disp('-----------------------------------------------------------')
% fprintf('\t\t\t STATE 1\t\tSTATE 2\n')
% fprintf('\t\t\t UNBURNED\t\tBURNED\n')
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

guess(1) = T2;
guess(2) = u1;
end

function strP = state(strR,r,T,phi,pP,E,S,C,M,PD,TN,strThProp)
% Calculate frozen state given T & rho
strR.v = strR.mi/r*1e3;
TP = T; % vP = vR (computed from R);
% Equilibrium composition at defined T and constant v
PD.ProblemType = 'TV';
strP = SolveProblemTP_TV(strR,phi,pP,TP,E,S,C,M,PD,TN,strThProp);
end