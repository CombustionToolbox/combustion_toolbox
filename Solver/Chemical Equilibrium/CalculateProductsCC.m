function [NCO2P_0,NCOP_0,NH2OP_0,NH2P_0,NO2P_0,NN2P_0,NHeP_0,NArP_0,NCgrP_0,phi_c,FLAG_SOOT] =  CalculateProductsCC(NatomE,phi,TP,Elements,factor_c,Fuel,strThProp)

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

% if x~= 0
%     phi_c = 2/(x)*(x+y/4);
% else
%     phi_c = phi;
% end
if Fuel.x~= 0 && (Fuel.x~= Fuel.z)
%     DG0 = (species_g0_new('CO2',TP,strThProp)-2*species_g0_new('CO',TP,strThProp))*1000;
%     k7 = exp(-DG0/(R0*TP));
%     Fuel.eps = 1/2*(-1 + sqrt(1 + 4*k7*Fuel.x));

    if Fuel.eps <= 0.5 && Fuel.eps>1e-16
%     phi_c = -(((-1+Fuel.eps)*(4*Fuel.x+Fuel.y-2*Fuel.z))/(2*(1+Fuel.eps)*Fuel.x+Fuel.eps^2*Fuel.y-2*Fuel.z+2*Fuel.eps*Fuel.z)); % C_x H_y O_z
    phi_c =(2*(Fuel.x+Fuel.y/4-Fuel.z/2))/(-Fuel.z+((2*Fuel.eps*Fuel.x+Fuel.x+Fuel.eps*Fuel.y/2)/(1+Fuel.eps))); % C_x H_y O_z
%     phi_c = -(((-1+Fuel.eps)*(4*Fuel.x+Fuel.y))/(2*(1+Fuel.eps)*Fuel.x+Fuel.eps^2*Fuel.y)); % C_x H_y
    else
        phi_c = 2/(Fuel.x-Fuel.z)*(Fuel.x+Fuel.y/4-Fuel.z/2);
    end
else
    phi_c = 1.1*phi;
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
%     elseif phi < phi_c
%     else
        % general case of rich mixtures with hydrogens (H) and carbons (C)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Equilibrium constant for the inverse wager-gas shift reaction
        %
        % CO2+H2 <-IV-> CO+H2O
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        DG0 = (species_g0_new('CO',TP,strThProp)+species_g0_new('H2O',TP,strThProp)-species_g0_new('CO2',TP,strThProp))*1000;
        k4 = exp(-DG0/(R0*TP));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        NCOP_0  = round((1/4)*(6*k4*x+k4*y-2*k4*z-4*x+2*z-sqrt(24*k4*x*z+16*x^2-16*x*z-16*k4*x^2+4*k4^2*x*y-8*k4^2*x*z-4*k4^2*y*z+4*k4*y*z+4*k4^2*x^2+k4^2*y^2+4*k4^2*z^2-8*k4*z^2+4*z^2))/(k4-1),14);
        
        NCO2P_0 =    x      -NCOP_0;
        NH2OP_0 = -2*x      +z+NCOP_0;
        NH2P_0  =  2*x+y/2-z-NCOP_0;
        
%     end
%     if  any(NCgrP_0 < 0 || NCO2P_0 < 0 || NCOP_0 < 0 || NH2OP_0 < 0 || NH2P_0 < 0 || NO2P_0 < 0 || NN2P_0 < 0 || NHeP_0 < 0 || NArP_0 < 0)
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
        
        DG0 = (species_g0_new('CO2',TP,strThProp)-2*species_g0_new('CO',TP,strThProp))*1000;
        k7 = exp(-DG0/(R0*TP));
        DG0 = (species_g0_new('CO',TP,strThProp)+species_g0_new('H2O',TP,strThProp)-species_g0_new('CO2',TP,strThProp))*1000;
        k4 = exp(-DG0/(R0*TP));
        

        zeta =1;
%         pP=1;
%         zeta = NP/pP;
        mu = k7/zeta;
%         NCOP_0 = real((1/(24*k4*mu^2))*(-4*(2+k4)*mu-(4*2^(1/3)*mu^2*(-4+2*k4+k4^2*(-1+3*y*mu-6*z*mu)))/(-16*mu^3+12*k4*mu^3+6*k4^2*mu^3-2*k4^3*mu^3+18*k4^2*y*mu^4+9*k4^3*y*mu^4+72*k4^2*z*mu^4-18*k4^3*z*mu^4+sqrt(mu^6*(4*(-4+2*k4+k4^2*(-1+3*y*mu-6*z*mu))^3+(-16+12*k4+k4^3*(-2+9*y*mu-18*z*mu)+6*k4^2*(1+3*y*mu+12*z*mu))^2)))^(1/3)+2*2^(2/3)*(-16*mu^3+12*k4*mu^3+6*k4^2*mu^3-2*k4^3*mu^3+18*k4^2*y*mu^4+9*k4^3*y*mu^4+72*k4^2*z*mu^4-18*k4^3*z*mu^4+sqrt(mu^6*(4*(-4+2*k4+k4^2*(-1+3*y*mu-6*z*mu))^3+(-16+12*k4+k4^3*(-2+9*y*mu-18*z*mu)+6*k4^2*(1+3*y*mu+12*z*mu))^2)))^(1/3)));
%         NCOP_0 = real(1/(48*k4*mu^2)*(-8*(2+k4)*mu+(4*2^(1/3)*(1+1i*sqrt(3))*mu^2*(-4+2*k4+k4^2*(-1+3*mu*(y-2*z))))/(-16*mu^3+12*k4*mu^3+6*k4^2*mu^3-2*k4^3*mu^3+18*k4^2*mu^4*y+9*k4^3*mu^4*y+72*k4^2*mu^4*z-18*k4^3*mu^4*z+sqrt(mu^6*(4*(-4+2*k4+k4^2*(-1+3*mu*(y-2*z)))^3+(-16+12*k4+k4^3*(-2+9*mu*(y-2*z))+6*k4^2*(1+3*mu*(y+4*z)))^2)))^(1/3)+2*1i*2^(2/3)*(1i+sqrt(3))*(-16*mu^3+12*k4*mu^3+6*k4^2*mu^3-2*k4^3*mu^3+18*k4^2*mu^4*y+9*k4^3*mu^4*y+72*k4^2*mu^4*z-18*k4^3*mu^4*z+sqrt(mu^6*(4*(-4+2*k4+k4^2*(-1+3*mu*(y-2*z)))^3+(-16+12*k4+k4^3*(-2+9*mu*(y-2*z))+6*k4^2*(1+3*mu*(y+4*z)))^2)))^(1/3)));
%         NCOP_0 = (sqrt(8*mu*z+1)-1)/(4*mu);
        

        a0 = -2*k4/zeta-k7*y+2*k7*z;
        a1 = 2*k7+4*k4*k7;
        a2 = 4*k7^2*zeta;
        a3 = 2*k4*z/zeta;
        
        NCOP_0 = real(-(a1/(3*a2))-(2^(1/3)*(-a1^2-3*a0*a2))/(3*a2*(-2*a1^3-9*a0*a1*a2+27*a2^2*a3+sqrt(-4*(a1^2+3*a0*a2)^3+(2*a1^3+9*a0*a1*a2-27*a2^2*a3)^2))^(1/3))+(-2*a1^3-9*a0*a1*a2+27*a2^2*a3+sqrt(-4*(a1^2+3*a0*a2)^3+(2*a1^3+9*a0*a1*a2-27*a2^2*a3)^2))^(1/3)/(3*2^(1/3)*a2));
%         a2 = log(4)+2*log(k7)+log(zeta);
%         NCOP_0 = real(-(a1/(3*exp(a2)))-(2^(1/3)*(-a1^2-3*a0*exp(a2)))/(3*exp(a2)*(-2*a1^3-9*a0*a1*exp(a2)+27*exp(a2)^2*a3+sqrt(-4*(a1^2+3*a0*exp(a2))^3+(2*a1^3+9*a0*a1*exp(a2)-27*exp(a2)^2*a3)^2))^(1/3))+(-2*a1^3-9*a0*a1*exp(a2)+27*exp(a2)^2*a3+sqrt(-4*(a1^2+3*a0*exp(a2))^3+(2*a1^3+9*a0*a1*exp(a2)-27*exp(a2)^2*a3)^2))^(1/3)/(3*2^(1/3)*exp(a2)));

        NCO2P_0 = mu*NCOP_0^2;
        NCgrP_0 = x-NCO2P_0-NCOP_0;
        NH2OP_0 = z-2*NCO2P_0-NCOP_0;
        NH2P_0  = y/2-NH2OP_0;
        
        FLAG_SOOT = 1;    
        
        %%%%%%%%%%%%%%%%%%%%%%% OTHER CASE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Equilibrium constant for the reaction
        %
        % CO2+2H2 <-VIII-> C(gr)+2H2O
        %
        % Equilibrium constant for the inverse wager-gas shift reaction
        %
        % CO2+H2 <-IV-> CO+H2O
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%         DG0 = (2*species_g0_new('H2O',TP,strThProp)-species_g0_new('CO2',TP,strThProp))*1000;
%         KPT_VIII = exp(-DG0/(R0*TP));
%         mu2 = zeta/KPT_VIII;
%         NH2OP_0 = (mu2+2*k4*mu2*(1+y)-sqrt(mu2*(mu2*(1+2*k4*(1+y))^2-16*k4^2*z)))/(8*k4*mu2);
%         NH2P_0  = y/2-NH2OP_0;
%         NCO2P_0 = mu2*(NH2OP_0/NH2P_0)^2;
%         NCOP_0 = z-2*NCO2P_0-NH2OP_0;
%         NCgrP_0 = x-NCO2P_0-NCOP_0;
        
    end
end