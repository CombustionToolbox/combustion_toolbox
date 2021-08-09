function phi_c = CalculatePhic(Fuel, Ninerts, phi, TP, pP, strThProp)
    R0 = 8.3144598; % [J/(K mol)]. Universal gas constant
    tol = 1;
    if Fuel.x~= 0 && (Fuel.x~= Fuel.z)
        DG0 = (species_g0('CO2',TP,strThProp)-2*species_g0('CO',TP,strThProp))*1000;
        k7 = exp(-DG0/(R0*TP));
        DG0 = (species_g0('CO',TP,strThProp)+species_g0('H2O',TP,strThProp)-species_g0('CO2',TP,strThProp))*1000;
        k4 = exp(-DG0/(R0*TP));
        
        phi_c0 = 2/(Fuel.x-Fuel.z)*(Fuel.x+Fuel.y/4-Fuel.z/2);
        phi_c = phi_c0;
        while tol > 1e-6 
            NP = Fuel.x+Fuel.y/2+Ninerts+Fuel.w/2+79/21/phi_c*(Fuel.x+Fuel.y/4-Fuel.z/2);
            zeta = NP/pP;
            nco = -((zeta - sqrt(zeta)*sqrt(4*k7*Fuel.x+zeta))/(2*k7));
            nco2 = Fuel.x-nco;
            nh2o = Fuel.y/(2*(1+nco/nco2*k4));
            
            phi_c_old = phi_c;
            phi_c = 2*(Fuel.x+Fuel.y/4-Fuel.z/2)/(2*nco2+nh2o+nco-Fuel.z);
            tol = abs((phi_c-phi_c_old)/phi_c);
        end
    else
        phi_c = 1.1*phi;
    end
end