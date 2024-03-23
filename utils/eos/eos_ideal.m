function V = eos_ideal(T, p, varargin)
    % Compute pressure considering ideal Equation of State (EoS)
    %
    % Args:
    %     T (float): Temperature of the mixture [K]
    %     p (float): Pressure of the mixture [Pa]
    % 
    % Returns:
    %     V (float): Molar volume of the mixture [m3/mol]
    
    R0 = 8.31446261815324; % [J/(K mol)]. Universal gas constant
    
    V = R0 * T ./ p;
end