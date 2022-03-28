function phi_c0 = Compute_phi_c(Fuel)
    % Compute guess of equivalence ratio in which soot appears considering complete combustion
    %
    % Args:
    %     Fuel (struct): Struct mix with all the properties of the Fuel mixture
    %
    % Returns:
    %     phi_c (float): Equivalence ratio in which soot appears [-]

    phi_c0 = 2/(Fuel.x - Fuel.z) * (Fuel.x + Fuel.y/4 - Fuel.z/2);
    if phi_c0 <= 1e-5 || isnan(phi_c0) 
       phi_c0 = inf; 
    end
end