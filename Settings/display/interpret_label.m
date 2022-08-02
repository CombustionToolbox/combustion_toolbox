function value = interpret_label(property)
    % Switch label from input property
    % 
    % Args:
    %     property (str): Property name
    %
    % Returns:
    %     value (str): Corresponding name of the property
    
    if iscell(property)
        if length(property) > 1
            value = 'Multiple variables';
            return
        end
        property = property{1};
    end

    switch lower(property)
        case 'phi'
            value = 'Equivalance ratio';
        case 'rho'
            value = 'Density [kg/m$^3$]';
        case 't'
            value = 'Temperature [K]';
        case 'p'
            value = 'Pressure [bar]';
        case 'h'
            value = 'Enthalpy [kJ/kg]';
        case 'e'
            value = 'Internal energy [kJ/kg]';
        case 'g'
            value = 'Gibbs energy [kJ/kg]';
        case {'s', 's0'}
            value = 'Entropy [kJ/kg-K]';
        case 'w'
            value = 'Molecular weight [g/mol]';
        case 'cp'
            value = '$c_p$ [kJ/kg-K]';
        case 'cv'
            value = '$c_v$ [kJ/kg-K]';
        case 'gamma_s'
            value = 'Adiabatic index';
        case {'sound', 'a'}
            value = 'Sound velocity [m/s]';
        case 'u'
             value = 'Incident velocity [m/s]';
        case 'v_shock'
            value = 'Shock velocity [m/s]';
        case 'u_preshock'
            value = 'Preshock velocity [m/s]';
        case 'u_postshock'
            value = 'Postshock velocity [m/s]';
        case 'dvdp_t'
            value = '$(\rm{d} v/ \rm{d} p)_T$';
        case 'dvdt_p'
            value = '$(\rm{d} v/ \rm{d} T)_p$';
        case 'cstar'
            value = 'Characteristic velocity [m/s]';
        case 'cf'
            value = 'Coeffficient of thrust';
        case 'i_sp'
            value = 'Specific impulse ambient [s]';
        case 'i_vac'
            value = 'Specific impulse vacuum [s]';
        case 'n'
            value = 'Total n$^\circ$ moles [mol]';
        case 'v'
            value = 'Volume [m$^3$]';
        case 'v_sp'
            value = 'Specific volume [m$^3$/kg]';
        case 'hf'
            value = 'Enthalpy of formation [kJ/kg]';
        case 'ef'
            value = 'Internal energy of formation [kJ/kg]';
        case 'error_moles'
            value = 'Relative error moles composition';
        case 'error_ions'
            value = 'Relative error electroneutrality';
        case 'error_problem'
            value = 'Relative error problem';
        otherwise
            value = ['$', property, '$'];
    end
end