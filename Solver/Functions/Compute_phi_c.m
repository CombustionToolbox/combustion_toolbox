function phi_c0 = Compute_phi_c(Fuel)
    phi_c0 = 2/(Fuel.x - Fuel.z) * (Fuel.x + Fuel.y/4 - Fuel.z/2);
end