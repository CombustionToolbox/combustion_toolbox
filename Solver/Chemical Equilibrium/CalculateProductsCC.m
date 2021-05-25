function [N_CC, phi_c, FLAG_SOOT] =  CalculateProductsCC(app, strR, phi, pP, TP)

% Abbreviations ---------------------
NatomE = strR.NatomE;
Elements = app.E.elements;
factor_c = app.TN.factor_c;
Fuel = app.PD.Fuel;
strThProp = app.strThProp;
% -----------------------------------

phi_c0 = Compute_phi_c(Fuel);

R0 = 8.3144598; % [J/(K mol)]. Universal gas constant

x = NatomE(strcmp(Elements,'C'));
y = NatomE(strcmp(Elements,'H'));
z = NatomE(strcmp(Elements,'O'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inerts
w = NatomE(strcmp(Elements,'N'));
b = NatomE(strcmp(Elements,'He'));
c = NatomE(strcmp(Elements,'Ar'));

NN2P_0 = w/2;
NHeP_0 = b;
NArP_0 = c;
NCgrP_0 = 0;

FLAG_SOOT = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ninerts = b + c;
phi_c = CalculatePhic(Fuel, Ninerts, phi, TP, pP, strThProp);
if phi_c0 <= 1e-5
   phi_c = 1e5; 
end
if phi <= 1 % case of lean or stoichiometric mixtures
    
    NCO2P_0 = x;
    NCOP_0  = 0;
    NH2OP_0 = y/2;
    NH2P_0  = 0;
    NO2P_0  =-x-y/4+z/2;
    
else % case of rich mixtures
    
    NO2P_0 = 0;
    
    if (x == 0) && (y ~= 0) % if there are only hydrogens (H)
        
        NCO2P_0 = 0;
        NCOP_0  = 0;
        NH2OP_0 = z;
        NH2P_0  = y/2-z;
        
    elseif (x ~= 0) && (y == 0) && phi < phi_c % if there are only carbons (C)
        
        NCO2P_0 = -x+z;
        NCOP_0  = 2*x-z;
        NH2OP_0 = 0;
        NH2P_0  = 0;
        
    elseif phi < phi_c*factor_c
        % general case of rich mixtures with hydrogens (H) and carbons (C)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Equilibrium constant for the inverse wager-gas shift reaction
        %
        % CO2+H2 <-IV-> CO+H2O
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        DG0 = (species_g0('CO',TP,strThProp)+species_g0('H2O',TP,strThProp)-species_g0('CO2',TP,strThProp))*1000;
        k4 = exp(-DG0/(R0*TP));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        NCOP_0  = round((1/4)*(6*k4*x+k4*y-2*k4*z-4*x+2*z-sqrt(24*k4*x*z+16*x^2-16*x*z-16*k4*x^2+4*k4^2*x*y-8*k4^2*x*z-4*k4^2*y*z+4*k4*y*z+4*k4^2*x^2+k4^2*y^2+4*k4^2*z^2-8*k4*z^2+4*z^2))/(k4-1),14);
        
        NCO2P_0 =    x      -NCOP_0;
        NH2OP_0 = -2*x      +z+NCOP_0;
        NH2P_0  =  2*x+y/2-z-NCOP_0;
        
    elseif phi >= phi_c*factor_c
        % general case of rich mixtures with hydrogens (H) carbons (C) and soot
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Equilibrium constant for the Boudouard reaction
        %
        % 2CO <-VII-> CO2+C(gr)
        %
        % Equilibrium constant for the inverse wager-gas shift reaction
        %
        % CO2+H2 <-IV-> CO+H2O
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        DG0 = (species_g0('CO2',TP,strThProp)-2*species_g0('CO',TP,strThProp))*1000;
        k7 = exp(-DG0/(R0*TP));
        DG0 = (species_g0('CO',TP,strThProp)+species_g0('H2O',TP,strThProp)-species_g0('CO2',TP,strThProp))*1000;
        k4 = exp(-DG0/(R0*TP));
        

        zeta =1;
        mu = k7/zeta;

        

        a0 = -2*k4/zeta-k7*y+2*k7*z;
        a1 = 2*k7+4*k4*k7;
        a2 = 4*k7^2*zeta;
        a3 = 2*k4*z/zeta;
        
        NCOP_0 = real(-(a1/(3*a2))-(2^(1/3)*(-a1^2-3*a0*a2))/(3*a2*(-2*a1^3-9*a0*a1*a2+27*a2^2*a3+sqrt(-4*(a1^2+3*a0*a2)^3+(2*a1^3+9*a0*a1*a2-27*a2^2*a3)^2))^(1/3))+(-2*a1^3-9*a0*a1*a2+27*a2^2*a3+sqrt(-4*(a1^2+3*a0*a2)^3+(2*a1^3+9*a0*a1*a2-27*a2^2*a3)^2))^(1/3)/(3*2^(1/3)*a2));
        
        NCO2P_0 = mu*NCOP_0^2;
        NCgrP_0 = x-NCO2P_0-NCOP_0;
        NH2OP_0 = z-2*NCO2P_0-NCOP_0;
        NH2P_0  = y/2-NH2OP_0;
        
        FLAG_SOOT = 1;    
    
    end
end
N_CC = [NCO2P_0, NCOP_0, NH2OP_0, NH2P_0, NO2P_0, NN2P_0, NHeP_0, NArP_0, NCgrP_0];