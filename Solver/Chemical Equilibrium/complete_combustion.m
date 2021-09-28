function [N0, species] = complete_combustion(self, str, phi)
% 
% Abbreviations ---------------------
Fuel = self.PD.Fuel;
% -----------------------------------
species = {'CO2', 'CO', 'H2O', 'H2', 'O2'};

if isempty(self.E.ind_C), x = 0; else, x = str.NatomE(self.E.ind_C); end
if isempty(self.E.ind_H), y = 0; else, y = str.NatomE(self.E.ind_H); end
if isempty(self.E.ind_O), z = 0; else, z = str.NatomE(self.E.ind_O); end

phi_c = Compute_phi_c(Fuel);

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
        disp('Not yet'); 
    elseif phi >= phi_c
        % general case of rich mixtures with hydrogens (H), carbons (C) and soot
        disp('Not yet');
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

function N0 = compute_moles_rich(x, y, z, T, R0, strThProp)
    DG0 = (species_g0('CO', T, strThProp) + species_g0('H2O', T, strThProp) - species_g0('CO2', T, strThProp))*1000;
    k4 = exp(-DG0 / (R0*T));

    NCO_0  = round((1/4)*(6*k4*x+k4*y-2*k4*z-4*x+2*z-sqrt(24*k4*x*z+16*x^2-16*x*z-16*k4*x^2+4*k4^2*x*y-8*k4^2*x*z-4*k4^2*y*z+4*k4*y*z+4*k4^2*x^2+k4^2*y^2+4*k4^2*z^2-8*k4*z^2+4*z^2))/(k4-1),14);
    NCO2_0 =    x - NCO_0;
    NH2O_0 = -2*x + z + NCO_0;
    NH2_0  =  2*x + y/2 - z - NCO_0;
    NO2_0  = 0;
    N0 = [NCO2_0, NCO_0, NH2O_0, NH2_0, NO2_0];
end
