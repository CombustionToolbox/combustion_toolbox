function p = eos_ideal_p(n, T, v, varargin)
    % Compute pressure considering ideal Equation of State (EoS)
    %
    % Args:
    %     n (float): Number of moles of the mixture in gaseous phase [mol]
    %     T (float): Temperature of the mixture [K]
    %     v (float): Volume of the mixture [m3]
    % 
    % Returns:
    %     p (float): Pressure of the mixture [Pa]
    
    R0 = 8.31446261815324; % [J/(K mol)]. Universal gas constant

    p = (n .* R0 .* T) ./ v;
end