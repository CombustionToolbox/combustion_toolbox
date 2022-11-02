function [funname_NASA, funname_CT, y_labelname] = set_inputs_thermo_validations(property)
    % Set corresponding thermodynamic functions for NASA and Combustion
    % Toolbox
    %
    % Args:
    %   property (str): Thermodynamic property name
    %
    % Returns:
    %     Tuple containing
    %
    %     - funname_NASA (function): Function to use NASA's polynomials
    %     - funname_CT (function): Function to use Combustion Toolbox polynomials
    %     - y_labelname (str): Label y axis

    switch lower(property)
        case 'cp'
            funname_NASA = @species_cP_NASA;
            funname_CT = @species_cP;
            y_labelname = 'Heat capacity at constant pressure [J/(mol-K)]';
        case 'cv'
            funname_NASA = @species_cV_NASA;
            funname_CT = @species_cV;
            y_labelname = 'Heat capacity at constant volume [J/(mol-K)]';
        case 'det'
            funname_NASA = @species_DeT_NASA;
            funname_CT = @species_DeT;
            y_labelname = 'Thermal internal energy [kJ/(mol)]';
        case 'dht'
            funname_NASA = @species_DhT_NASA;
            funname_CT = @species_DhT;
            y_labelname = 'Thermal enthalpy [kJ/(mol)]';
        case 'g0'
            funname_NASA = @species_g0_NASA;
            funname_CT = @species_g0;
            y_labelname = 'Gibbs energy [kJ/(mol)]';
        case 'h0'
            funname_NASA = @species_h0_NASA;
            funname_CT = @species_h0;
            y_labelname = 'Enthalpy [kJ/(mol)]';
        case 's0'
            funname_NASA = @species_s0_NASA;
            funname_CT = @species_s0;
            y_labelname = 'Entropy [kJ/(mol-K)]';
        otherwise
            error('There is not such property on files');
    end

end
