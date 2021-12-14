function phi_c0 = Compute_phi_c(Fuel)
    phi_c0 = 2/(Fuel.x - Fuel.z) * (Fuel.x + Fuel.y/4 - Fuel.z/2);
    if phi_c0 <= 1e-5 || isnan(phi_c0) 
       phi_c0 = inf; 
    end
end