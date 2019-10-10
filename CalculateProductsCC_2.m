function [NCO2P_0,NCOP_0,NH2OP_0,NH2P_0,NO2P_0,NN2P_0,NHeP_0,NArP_0,NCgrP_0,phi_c,FLAG_SOOT] =  CalculateProductsCC_2(NatomE,phi,TP,Elements,factor_c,strThProp)

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
% phi_c = 2/(x-z)*(x+y/4-z/2);
if x~= 0
    phi_c = 2/(x)*(x+y/4);
else
    phi_c = phi;
end
if phi <= 1 % case of lean or stoichiometric mixtures
    
    NCO2P_0 = x;
    NCOP_0  = 0;
    NH2OP_0 =       y/2;
    NH2P_0  = 0;
    NO2P_0  =-x-y/4+z/2;
    
else % case of rich mixtures
    
    NO2P_0 = 0;
    
    if (x == 0) && (y ~= 0) % if there are only hydrogens (H)
        
        NCO2P_0 = 0;
        NCOP_0  = 0;
        NH2OP_0 =       z;
        NH2P_0  = y/2-z;
        
    elseif (x ~= 0) && (y == 0) && phi < phi_c % if there are only carbons (C)
        
        NCO2P_0 =-x+z;
        NCOP_0  = 2*x-z;
        NH2OP_0 = 0;
        NH2P_0  = 0;
        
%     elseif phi < phi_c*factor_c
%     elseif phi < phi_c
    else
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
        
        DG0 = (species_g0_new('CO2',TP,strThProp)-2*species_g0_new('CO',TP,strThProp))*1000;
        KPT_VII = exp(-DG0/(R0*TP));
        DG0 = (species_g0_new('CO',TP,strThProp)+species_g0_new('H2O',TP,strThProp)-species_g0_new('CO2',TP,strThProp))*1000;
        KPT_IV = exp(-DG0/(R0*TP));
        
        a0 = -2*KPT_IV-KPT_VII*y+2*KPT_VII*z;
        a1 = 2*KPT_VII+4*KPT_IV*KPT_VII;
        a2 = 4*KPT_VII^2;
        a3 = 2*KPT_IV*z;
        
        NCOP_0 = real(-(a1/(3*a2))-(2^(1/3)*(-a1^2-3*a0*a2))/(3*a2*(-2*a1^3-9*a0*a1*a2+27*a2^2*a3+sqrt(-4*(a1^2+3*a0*a2)^3+(2*a1^3+9*a0*a1*a2-27*a2^2*a3)^2))^(1/3))+(-2*a1^3-9*a0*a1*a2+27*a2^2*a3+sqrt(-4*(a1^2+3*a0*a2)^3+(2*a1^3+9*a0*a1*a2-27*a2^2*a3)^2))^(1/3)/(3*2^(1/3)*a2));
        
        NCO2P_0 = KPT_VII*NCOP_0^2;
        NCgrP_0 = x-NCO2P_0-NCOP_0;
        NH2OP_0 = z-2*NCO2P_0-NCOP_0;
        NH2P_0  = y/2-NH2OP_0;
        
        FLAG_SOOT = 1;
    end
    if  any(NCgrP_0 < 0 || NCO2P_0 < 0 || NCOP_0 < 0 || NH2OP_0 < 0 || NH2P_0 < 0 || NO2P_0 < 0 || NN2P_0 < 0 || NHeP_0 < 0 || NArP_0 < 0)
        % general case of rich mixtures with hydrogens (H) and carbons (C)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Equilibrium constant for the inverse wager-gas shift reaction
        %
        % CO2+H2 <-IV-> CO+H2O
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        DG0 = (species_g0_new('CO',TP,strThProp)+species_g0_new('H2O',TP,strThProp)-species_g0_new('CO2',TP,strThProp))*1000;
        KPT_IV = exp(-DG0/(R0*TP));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        NCOP_0  = round((1/4)*(6*KPT_IV*x+KPT_IV*y-2*KPT_IV*z-4*x+2*z-sqrt(24*KPT_IV*x*z+16*x^2-16*x*z-16*KPT_IV*x^2+4*KPT_IV^2*x*y-8*KPT_IV^2*x*z-4*KPT_IV^2*y*z+4*KPT_IV*y*z+4*KPT_IV^2*x^2+KPT_IV^2*y^2+4*KPT_IV^2*z^2-8*KPT_IV*z^2+4*z^2))/(KPT_IV-1),14);
        
        NCO2P_0 =    x      -NCOP_0;
        NH2OP_0 = -2*x      +z+NCOP_0;
        NH2P_0  =  2*x+y/2-z-NCOP_0;
        NCgrP_0 = 0;
        
        FLAG_SOOT = 0;
    end
end