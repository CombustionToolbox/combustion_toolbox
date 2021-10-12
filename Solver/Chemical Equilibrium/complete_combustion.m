function [N0, species] = complete_combustion(self, mix, phi)
% 
% Parameters ---------------------
Fuel = self.PD.Fuel;
phi_c = Compute_phi_c(Fuel);
% -----------------------------------
species = {'CO2', 'CO', 'H2O', 'H2', 'O2'};

if isempty(self.E.ind_C), x = 0; phi_c = inf; else, x = mix.NatomE(self.E.ind_C); end
if isempty(self.E.ind_H), y = 0; else, y = mix.NatomE(self.E.ind_H); end
if isempty(self.E.ind_O), z = 0; else, z = mix.NatomE(self.E.ind_O); end

if phi <= 1
    % case of lean or stoichiometric mixtures
    N0 = compute_moles_1ean(x, y, z);    
else
    % case of rich mixtures
    if (x == 0) && (y ~= 0)
        % if there are only hydrogens (H)
        N0 = compute_moles_rich_hydrogen(y, z); 
    elseif (x ~= 0) && (y == 0) && phi < phi_c
        % if there are only carbons (C)
        N0 = compute_moles_rich_carbon(x, z); 
    elseif phi < phi_c
        % general case of rich mixtures with hydrogens (H) and carbons (C)
%         N0 = compute_moles_rich_appr(x, y, z);
        T  = 1500;
        N0 = compute_moles_rich(x, y, z, T, self.C.R0, self.DB);
    elseif phi >= phi_c
        % general case of rich mixtures with hydrogens (H), carbons (C) and soot
        T  = 1000;
        N0 = compute_moles_rich_soot(y, z, T, self.C.R0, self.DB);
    end
end
end

% NESTED FUNCTIONS
function N0 = compute_moles_1ean(x, y ,z)
    NCO2_0 = x;
    NCO_0  = 0;
    NH2O_0 = y/2;
    NH2_0  = 0;
    NO2_0  = -x - y/4 + z/2;
    N0 = [NCO2_0, NCO_0, NH2O_0, NH2_0, NO2_0];
end

function N0 = compute_moles_rich_hydrogen(y, z)
    NCO2_0 = 0;
    NCO_0  = 0;
    NH2O_0 = z;
    NH2_0  = y/2 - z;
    NO2_0  = 0;
    N0 = [NCO2_0, NCO_0, NH2O_0, NH2_0, NO2_0];
end

function N0 = compute_moles_rich_carbon(x, z)
    NCO2_0 = -x + z;
    NCO_0  = 2*x - z;
    NH2O_0 = 0;
    NH2_0  = 0;
    NO2_0  = 0;
    N0 = [NCO2_0, NCO_0, NH2O_0, NH2_0, NO2_0];
end

function N0 = compute_moles_rich_appr(x, y, z)
    z = z/2;
    NCO2_0 = -x + z;
    NCO_0  = 2*x - z;
    NH2O_0 = z;
    NH2_0  = y/2 - z;
    NO2_0  = 0;
    N0 = [NCO2_0, NCO_0, NH2O_0, NH2_0, NO2_0];
end

function N0 = compute_moles_rich(x, y, z, T, R0, strThProp)
    DG0 = (species_g0('CO', T, strThProp) + species_g0('H2O', T, strThProp) - species_g0('CO2', T, strThProp)) * 1000;
    k4  = exp(-DG0 / (R0*T));

    NCO_0  = round((1/4)*(6*k4*x+k4*y-2*k4*z-4*x+2*z-sqrt(24*k4*x*z+16*x^2-16*x*z-16*k4*x^2+4*k4^2*x*y-8*k4^2*x*z-4*k4^2*y*z+4*k4*y*z+4*k4^2*x^2+k4^2*y^2+4*k4^2*z^2-8*k4*z^2+4*z^2))/(k4-1),14);
    NCO2_0 =    x - NCO_0;
    NH2O_0 = -2*x + z + NCO_0;
    NH2_0  =  2*x + y/2 - z - NCO_0;
    NO2_0  = 0;
    N0 = [NCO2_0, NCO_0, NH2O_0, NH2_0, NO2_0];
end

function N0 = compute_moles_rich_soot(y, z, T, R0, strThProp)
    DG0 = (species_g0('CO2', T, strThProp) - 2*species_g0('CO', T, strThProp)) * 1000;
    k7  = exp(-DG0 / (R0*T));
    DG0 = (species_g0('CO', T, strThProp) + species_g0('H2O', T, strThProp) - species_g0('CO2', T, strThProp)) * 1000;
    k4  = exp(-DG0 / (R0*T));

    zeta =1;
    mu = k7/zeta;

    a0 = -2*k4 / zeta - k7 * y + 2*k7 * z;
    a1 = 2*k7 + 4*k4 * k7;
    a2 = 4*k7^2 * zeta;
    a3 = 2*k4 * z / zeta;

    NCO_0 = real(-(a1/(3*a2))-(2^(1/3)*(-a1^2-3*a0*a2))/(3*a2*(-2*a1^3-9*a0*a1*a2+27*a2^2*a3+sqrt(-4*(a1^2+3*a0*a2)^3+(2*a1^3+9*a0*a1*a2-27*a2^2*a3)^2))^(1/3))+(-2*a1^3-9*a0*a1*a2+27*a2^2*a3+sqrt(-4*(a1^2+3*a0*a2)^3+(2*a1^3+9*a0*a1*a2-27*a2^2*a3)^2))^(1/3)/(3*2^(1/3)*a2));

    NCO2_0 = mu * NCO_0^2;
%     NCgr_0 = x - NCO2_0 - NCO_0; % Not required for our computations
    NH2O_0 = z - 2*NCO2_0 - NCO_0;
    NH2_0  = y/2 - NH2O_0;
    NO2_0  = 0;
    
    N0 = [NCO2_0, NCO_0, NH2O_0, NH2_0, NO2_0];
end
