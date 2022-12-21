function Gammas = compute_Gammas(u2, rho2, p2)
    % Compute slope of the Hugoniot curve
    %
    % Args:
    %     u2 (float): Post-shock velocity [m/s]
    %     rho2 (float): Post-shock density [kg/m3]
    %     p2 (float): Post-shock pressure [bar]
    %
    % Returns:
    %     Gammas (float): Slope of the Hugoniot curve [-]
    
    p2 = convert_bar_to_Pa(p2);
    
    Gammas =  u2(1:end-1).^2 .* compute_first_derivative(rho2, p2);
end