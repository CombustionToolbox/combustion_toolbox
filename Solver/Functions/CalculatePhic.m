function phi_c = CalculatePhic(Fuel, Ninerts, phi, TP, pP, DB)
    % Compute equivalence ratio in which soot appears considering complete combustion
    %
    % Args:
    %     Fuel (struct):     Struct mix with all the properties of the Fuel mixture
    %     Ninerts (float):   Number of moles of the inerts species
    %     phi (float):       Equivalence ratio [-]
    %     TP (float):        Temperature [K]
    %     pP (float):        Pressure [bar]
    %     DB (struct):       Database
    %
    % Returns:
    %     phi_c (float):     Equivalence ratio in which soot appears [-]

    R0 = 8.3144598; % [J/(K mol)]. Universal gas constant
    tol = 1;
    if Fuel.x~= 0 && (Fuel.x~= Fuel.z)
        DG0 = (species_g0('CO2', TP, DB) - 2*species_g0('CO', TP, DB))*1000;
        k7 = exp(-DG0 / (R0*TP));
        DG0 = (species_g0('CO', TP, DB) + species_g0('H2O', TP, DB) - species_g0('CO2', TP, DB))*1000;
        k4 = exp(-DG0 / (R0*TP));
        
        phi_c0 = 2/(Fuel.x-Fuel.z)*(Fuel.x+Fuel.y/4-Fuel.z/2);
        phi_c = phi_c0;
        while tol > 1e-6 
            NP = Fuel.x + Fuel.y/2 + Ninerts + Fuel.w/2 + 79/21/phi_c*(Fuel.x + Fuel.y/4 - Fuel.z/2);
            zeta = NP/pP;
            nco = -((zeta - sqrt(zeta) * sqrt(4*k7*Fuel.x + zeta)) / (2*k7));
            nco2 = Fuel.x - nco;
            nh2o = Fuel.y / (2*(1+nco/nco2*k4));
            
            phi_c_old = phi_c;
            phi_c = 2*(Fuel.x + Fuel.y/4 - Fuel.z/2) / (2*nco2+nh2o+nco-Fuel.z);
            tol = abs((phi_c - phi_c_old) / phi_c);
        end
    else
        phi_c = 1.1*phi;
    end
end