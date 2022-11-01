function value = interpreter_label(property, varargin)
    % Interpreter label for properties - returns property name
    %
    % Args:
    %     property (str): Property name
    %
    % Returns:
    %     value (str): Corresponding name of the property

    % Default values
    type = 'medium';

    
    % Check if property is a cell variable
    if iscell(property)

        if length(property) > 1
            value = 'Multiple variables';
            return
        end

        property = property{1};
    end
    
    % Check if property is empty
    if isempty(property)
        value = '';
        return
    end
    
    % Check additional inputs
    if nargin > 1
        type = varargin{1};
    end

    % Get labels
    [property_name, property_latex, property_unit] = property_names(property, type);

    % Convert property definition to its name
    switch lower(type)
        case 'short'

            if ~isempty(property_latex)
                property_name = '';
            end

        case 'medium'

            if ~isempty(property_name)
                property_latex = '';
            end

        case 'long'
            % nothing
    end

    value = sprintf('%s %s %s', property_name, property_latex, property_unit);
end

% SUB-PASS FUNCTIONS
function [property_name, property_latex, property_unit] = property_names(property, type)

    switch lower(property)
        case 'phi'
            property_name = 'Equivalance ratio';
            property_latex = '\phi';
            property_unit = '';
        case 'rho'
            property_name = 'Density';
            property_latex = '\rho';
            property_unit = '[kg/m$^3$]';
        case 't'
            property_name = 'Temperature';
            property_latex = 'T';
            property_unit = '[K]';
        case 'p'
            property_name = 'Pressure';
            property_latex = 'p';
            property_unit = '[bar]';
        case 'h'
            property_name = 'Enthalpy';
            property_latex = 'h';
            property_unit = '[kJ/kg]';
        case 'e'
            property_name = 'Internal energy';
            property_latex = 'e';
            property_unit = '[kJ/kg]';
        case 'g'
            property_name = 'Gibbs energy';
            property_latex = 'g';
            property_unit = '[kJ/kg]';
        case 's'
            property_name = 'Entropy';
            property_latex = 's';
            property_unit = '[kJ/kg-K]';
        case 's0'
            property_name = 'Entropy frozen';
            property_latex = 's_0';
            property_unit = '[kJ/kg-K]';
        case 'ds'
            property_name = 'Entropy of mixing';
            property_latex = '\Delta s';
            property_unit = '[kJ/kg-K]';
        case 'w'
            property_name = 'Molecular weight';
            property_latex = 'W';
            property_unit = '[g/mol]';
        case 'cp'
            property_name = 'Specific heat pressure';
            property_latex = 'c_p';
            property_unit = '[kJ/kg-K]';
        case 'cp_f'
            property_name = 'Specific heat pressure frozen';
            property_latex = 'c_{p,f}';
            property_unit = '[kJ/kg-K]';
        case 'cp_r'
            property_name = 'Specific heat pressure reaction';
            property_latex = 'c_{p,r}';
            property_unit = '[kJ/kg-K]';
        case 'cv'
            property_name = 'Specific heat volume';
            property_latex = 'c_v';
            property_unit = '[kJ/kg-K]';
        case 'cv_f'
            property_name = 'Specific heat volume frozen';
            property_latex = 'c_{v,f}';
            property_unit = '[kJ/kg-K]';
        case 'cv_r'
            property_name = 'Specific heat volume reaction';
            property_latex = 'c_{v,r}';
            property_unit = '[kJ/kg-K]';
        case 'gamma'
            property_name = 'Specific heat ratio';
            property_latex = '\gamma';
            property_unit = '';
        case 'gamma_s'
            property_name = 'Adibatic index';
            property_latex = '\gamma_s';
            property_unit = '';
        case {'sound', 'a'}
            property_name = 'Sound velocity';
            property_latex = 'a';
            property_unit = '[m/s]';
        case 'u'
            property_name = 'Velocity';
            property_latex = 'u';
            property_unit = '[m/s]';
        case 'v_shock'
            property_name = 'Shock velocity';
            property_latex = 'v_{\rm shock}';
            property_unit = '[m/s]';
        case 'u_preshock'
            property_name = 'Pre-shock velocity';
            property_latex = 'u_{\rm preshock}';
            property_unit = '[m/s]';
        case 'u_postshock'
            property_name = 'Post-shock velocity';
            property_latex = 'u_{\rm postshock}';
            property_unit = '[m/s]';
        case 'cstar'
            property_name = 'Characteristic velocity';
            property_latex = 'C^*';
            property_unit = '[m/s]';
        case 'cf'
            property_name = 'Coefficient of thrust';
            property_latex = 'C_f';
            property_unit = '';
        case 'i_sp'
            property_name = 'Specific impulse ambient';
            property_latex = 'I_{sp}';
            property_unit = '[s]';
        case 'i_vac'
            property_name = 'Specific impulse vaccum';
            property_latex = 'I_{vac}';
            property_unit = '[s]';
        case 'n'
            property_name = 'Total moles';
            property_latex = 'n';
            property_unit = '[mol]';
        case 'v'
            property_name = 'Volume';
            property_latex = 'v';
            property_unit = '[m$^3$]';
        case 'v_sp'
            property_name = 'Specific volume';
            property_latex = 'v_{sp}';
            property_unit = '[m$^3$/kg]';
        case 'hf'
            property_name = 'Enthalpy of formation';
            property_latex = 'h_f';
            property_unit = '[kJ/kg]';
        case 'dht'
            property_name = 'Enthalpy thermal';
            property_latex = '\Delta h_T';
            property_unit = '[kJ/kg]';
        case 'ef'
            property_name = 'Internal energy of formation';
            property_latex = 'e_f';
            property_unit = '[kJ/kg]';
        case 'det'
            property_name = 'Internal energy thermal';
            property_latex = '\Delta e_T';
            property_unit = '[kJ/kg]';
        case 'xi'
            property_name = 'Molar fractions';
            property_latex = 'X_i';
            property_unit = '';
        case 'yi'
            property_name = 'Mass fractions';
            property_latex = 'Y_i';
            property_unit = '';
        case 'v_p/v_r'
            property_name = 'Volume ratio';
            property_latex = 'v_2/v_1';
            property_unit = '';
        case 'error_moles'
            property_name = 'Relative error moles composition';
            property_latex = '\epsilon_{\rm moles}';
            property_unit = '';
        case 'error_moles_ions'
            property_name = 'Relative error electroneutrality';
            property_latex = '\epsilon_{\rm ions}';
            property_unit = '';
        case 'error_problem'
            property_name = 'Relative error problem';
            property_latex = '\epsilon_{\rm problem}';
            property_unit = '';
        case {'dvdtp', 'dvdt_p'}
            property_name = '';
            property_latex = '(\rm{d}V/\rm{d}T)_p';
            property_unit = '';
        case {'dvdpt', 'dvdp_t'}
            property_name = '';
            property_latex = '(\rm{d}V/\rm{d}p)_T';
            property_unit = '';
        case 'theta'
            property_name = 'Deflection angle';
            property_latex = '\theta';
            property_unit = '[deg]';
        case 'beta'
            property_name = 'Wave angle';
            property_latex = '\beta';
            property_unit = '[deg]';
        case 'overdriven'
            property_name = 'Overdriven ratio';
            property_latex = 'u_1/u_{\rm cj}';
            property_unit = '';
        case 'underdriven'
            property_name = 'Underdriven ratio';
            property_latex = 'u_1/u_{\rm cj}';
            property_unit = '';
        case 'drive_factor'
            property_name = 'Drive factor';
            property_latex = 'u_1/u_{\rm cj}';
            property_unit = '';
        case 'of'
            property_name = 'Mixture ratio';
            property_latex = 'O/F';
            property_unit = '';
        case 'mi'
            property_name = 'Mass';
            property_latex = 'm';
            property_unit = '[kg]';
        case 'pv'
            property_name = 'Pressure $\times$ Volume';
            property_latex = 'pv';
            property_unit = '[bar-m$^3$]';
        otherwise
            property_unit = '';

            switch type
                case 'short'
                    property_name = '';
                    property_latex = property;
                case {'medium', 'long'}
                    property_name = property;
                    property_latex = '';
            end

    end

    property_latex = ['$', property_latex, '$'];

    if strcmpi(type, 'long')
        property_name = [property_name, ','];
    end

end
