function value = interpreterLabel(property, varargin)
    % Interpreter label for properties - returns property name
    %
    % Note:
    %     The 'interpreter_label.m' routine considers that the properties
    %     are in mass basis. This will be fixed in a future patch.
    %
    % Args:
    %     property (char): Property name
    %
    % Optional Args:
    %     * type (char): Type of label to return. Can be 'short', 'medium' or 'long' (default: medium)
    %     * FLAG_BASIS (bool): Flag indicating label with basis, e.g., kg or mol (default: true)
    %     * basis (char): Specific basis (default: kg)
    %
    % Returns:
    %     value (char): Corresponding name of the property

    % Default values
    type = 'medium';
    FLAG_BASIS = true;
    basis = 'kg';

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
    for i = 1:nargin-1

        switch i
            case 1
                type = varargin{1};
            case 2
                FLAG_BASIS = varargin{2};
            case 3
                basis = check_basis(varargin{3});
        end

    end

    

    % Get labels
    [property_name, property_latex, property_unit] = property_names(property, type, FLAG_BASIS, basis);

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
function basis = check_basis(basis)
    % Check basis

    switch lower(basis)
        case 'mi'
            basis = 'kg';
        case 'mw'
            basis = 'mol';
        otherwise
            error('Not known basis.');
    end

end

function [property_name, property_latex, property_unit] = property_names(property, type, FLAG_BASIS, basis)

    switch lower(property)
        case {'phi', 'equivalenceratio'}
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

            if FLAG_BASIS
                property_unit = sprintf('[kJ/%s]', basis);
            else
                property_unit = '[kJ]';
            end

        case 'e'
            property_name = 'Internal energy';
            property_latex = 'e';

            if FLAG_BASIS
                property_unit = sprintf('[kJ/%s]', basis);
            else
                property_unit = '[kJ]';
            end

        case 'g'
            property_name = 'Gibbs energy';
            property_latex = 'g';
            
            if FLAG_BASIS
                property_unit = sprintf('[kJ/%s]', basis);
            else
                property_unit = '[kJ]';
            end

        case 's'
            property_name = 'Entropy';
            property_latex = 's';
            
            if FLAG_BASIS
                property_unit = sprintf('[kJ/%s-K]', basis);
            else
                property_unit = '[kJ/K]';
            end

        case 's0'
            property_name = 'Entropy frozen';
            property_latex = 's_0';
            
            if FLAG_BASIS
                property_unit = sprintf('[kJ/%s-K]', basis);
            else
                property_unit = '[kJ/K]';
            end

        case 'ds'
            property_name = 'Entropy of mixing';
            property_latex = '\Delta s';
            
            if FLAG_BASIS
                property_unit = sprintf('[kJ/%s-K]', basis);
            else
                property_unit = '[kJ/K]';
            end

        case 'w'
            property_name = 'Molecular weight';
            property_latex = 'W';
            property_unit = '[g/mol]';
        case 'cp'
            property_name = 'Specific heat pressure';
            property_latex = 'c_p';

            if FLAG_BASIS
                property_unit = sprintf('[kJ/%s-K]', basis);
            else
                property_unit = '[kJ/K]';
            end

        case 'cp_f'
            property_name = 'Specific heat pressure frozen';
            property_latex = 'c_{p,f}';

            if FLAG_BASIS
                property_unit = sprintf('[kJ/%s-K]', basis);
            else
                property_unit = '[kJ/K]';
            end

        case 'cp_r'
            property_name = 'Specific heat pressure reaction';
            property_latex = 'c_{p,r}';

            if FLAG_BASIS
                property_unit = sprintf('[kJ/%s-K]', basis);
            else
                property_unit = '[kJ/K]';
            end

        case 'cv'
            property_name = 'Specific heat volume';
            property_latex = 'c_v';

            if FLAG_BASIS
                property_unit = sprintf('[kJ/%s-K]', basis);
            else
                property_unit = '[kJ/K]';
            end

        case 'cv_f'
            property_name = 'Specific heat volume frozen';
            property_latex = 'c_{v,f}';

            if FLAG_BASIS
                property_unit = sprintf('[kJ/%s-K]', basis);
            else
                property_unit = '[kJ/K]';
            end

        case 'cv_r'
            property_name = 'Specific heat volume reaction';
            property_latex = 'c_{v,r}';

            if FLAG_BASIS
                property_unit = sprintf('[kJ/%s-K]', basis);
            else
                property_unit = '[kJ/K]';
            end
            
        case 'gamma'
            property_name = 'Specific heat ratio';
            property_latex = '\gamma';
            property_unit = '';
        case 'gamma_f'
            property_name = 'Frozen specific heat ratio';
            property_latex = '\gamma_f';
            property_unit = '';
        case 'gamma_s'
            property_name = 'Adiabatic index';
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
        case {'v_shock', 'ushock'}
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
        case {'m', 'mach'}
            property_name = 'Mach number';
            property_latex = '\mathcal{M}';
            property_unit = '';
        case 'm1'
            property_name = 'Pre-shock Mach number';
            property_latex = '\mathcal{M}_1';
            property_unit = '';
        case 'm2'
            property_name = 'Post-shock Mach number';
            property_latex = '\mathcal{M}_2';
            property_unit = '';
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
        case {'v_sp', 'vspecific'}
            property_name = 'Specific volume';
            property_latex = 'v_{sp}';
            property_unit = '[m$^3$/kg]';
        case 'hf'
            property_name = 'Enthalpy of formation';
            property_latex = 'h_f';
            
            if FLAG_BASIS
                property_unit = sprintf('[kJ/%s]', basis);
            else
                property_unit = '[kJ]';
            end

        case 'dht'
            property_name = 'Enthalpy thermal';
            property_latex = '\Delta h_T';
            
            if FLAG_BASIS
                property_unit = sprintf('[kJ/%s]', basis);
            else
                property_unit = '[kJ]';
            end

        case 'ef'
            property_name = 'Internal energy of formation';
            property_latex = 'e_f';
            
            if FLAG_BASIS
                property_unit = sprintf('[kJ/%s]', basis);
            else
                property_unit = '[kJ]';
            end

        case 'det'
            property_name = 'Internal energy thermal';
            property_latex = '\Delta e_T';
            
            if FLAG_BASIS
                property_unit = sprintf('[kJ/%s]', basis);
            else
                property_unit = '[kJ]';
            end

        case 'xi'
            property_name = 'Molar fractions';
            property_latex = 'X_j';
            property_unit = '';
        case 'yi'
            property_name = 'Mass fractions';
            property_latex = 'Y_j';
            property_unit = '';
        case {'v_p/v_r', 'v_v', 'vp_vr'}
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
            property_latex = '(\rm{d}v/\rm{d}T)_p';
            property_unit = '';
        case {'dvdpt', 'dvdp_t'}
            property_name = '';
            property_latex = '(\rm{d}v/\rm{d}p)_T';
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
        case {'drive_factor', 'drivefactor'}
            property_name = 'Overdriven factor';
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
        case {'aratio', 'arearatio'}
            property_name = 'Area exit / throat';
            property_latex = 'A_{\rm ratio} = A_e/A_t';
            property_unit = '';
        case {'aratio_c', 'arearatiochamber'}
            property_name = 'Area combustor / throat';
            property_latex = 'A_{\rm ratio, c} = A_c/A_t';
            property_unit = '';
        case {'density_variance', 'rho_variance'}
            property_name = 'Density variance';
            property_latex = '\rho_{\rm rms}^2 / \rho_1^2';
            property_unit = '';
        case {'pressure_variance', 'p_variance'}
            property_name = 'Pressure variance';
            property_latex = 'p_{\rm rms}^2 / (\rho_1 a_1^2)^2';
            property_unit = '';
        case {'temperature_variance', 't_variance'}
            property_name = 'Temperature variance';
            property_latex = 'T_{\rm rms}^2 / T_1^2';
            property_unit = '';
        case {'velocity_variance', 'u_variance'}
            property_name = 'Velocity variance';
            property_latex = 'u_{\rm rms}^2 / a_1^2';
            property_unit = '';
        case {'v_variance'}
            property_name = 'Velocity variance';
            property_latex = 'v_{\rm rms}^2 / a_1^2';
            property_unit = '';
        case {'w_variance'}
            property_name = 'Velocity variance';
            property_latex = 'w_{\rm rms}^2 / a_1^2';
            property_unit = '';
        case {'corr_rhot'}
            property_name = 'Correlation density-temperature';
            property_latex = '\langle \rho'' T'' \rangle / (\rho_1 T_1)';
            property_unit = '';
        case {'corr_rhop'}
            property_name = 'Correlation density-pressure';
            property_latex = '\langle \rho'' p'' \rangle / (\rho_1 p_1)';
            property_unit = '';
        case {'corr_rhou'}
            property_name = 'Correlation density-velocity';
            property_latex = '\langle \rho'' u'' \rangle / (\rho_1 a_1)';
            property_unit = '';
        case {'corr_rhov'}
            property_name = 'Correlation density-velocity';
            property_latex = '\langle \rho'' v'' \rangle / (\rho_1 a_1)';
            property_unit = '';
        case {'corr_rhow'}
            property_name = 'Correlation density-velocity';
            property_latex = '\langle \rho'' w'' \rangle / (\rho_1 a_1)';
            property_unit = '';
        case {'r11', 'r11_f', 'r11_tke', 'r11_ftke'}
            property_name = 'Streamwise Reynolds stress';
            property_latex = 'R_{11}';
            property_unit = '';
        case {'r22', 'r22_f', 'r22_tke', 'r22_ftke'}
            property_name = 'Transverse Reynolds stress';
            property_latex = 'R_{22}';
            property_unit = '';
        case {'r33', 'r33_f', 'r33_tke', 'r33_ftke'}
            property_name = 'Transverse Reynolds stress';
            property_latex = 'R_{33}';
            property_unit = '';
        case {'rtt', 'rtt_f', 'rtt_tke', 'rtt_ftke'}
            property_name = 'Transverse Reynolds stress';
            property_latex = 'R_{\rm TT}';
            property_unit = '';
        case {'k', 'tke', 'tke_f'}
            property_name = 'Turbulent kinetic energy';
            property_latex = 'K';
            property_unit = '';
        case {'enstrophy33'}
            property_name = 'z-component of enstrophy';
            property_latex = 'W_z';
            property_unit = '';
        case {'enstrophytt'}
            property_name = 'Transverse enstrophy';
            property_latex = 'W_\perp';
            property_unit = '';
        case {'enstrophy'}
            property_name = 'Enstrophy';
            property_latex = 'W';
            property_unit = '';
        case {'lengthkolmogorov'}
            property_name = 'Kolmogorov length';
            property_latex = '\ell_{k}';
            property_unit = '';
        case {'ratiolengthkolmogorov'}
            property_name = 'Normalized Kolmogorov length';
            property_latex = '\ell_{k} / \ell_{k, 1}';
            property_unit = '';
        case {'ratiogridlengthkolmogorov'}
            property_name = 'Ratio grid spacing to Kolmogorov length';
            property_latex = '\Delta x / \ell_{k}';
            property_unit = '';
        case {'rratio'}
            property_name = 'Density ratio';
            property_latex = '\mathcal{R}';
            property_unit = '';
        case {'pratio'}
            property_name = 'Pressure ratio';
            property_latex = '\mathcal{P}';
            property_unit = '';
        case {'tratio'}
            property_name = 'Temperature ratio';
            property_latex = '\mathcal{T}';
            property_unit = '';
        case {'gammas1'}
            property_name = 'Dimensionless slope RH';
            property_latex = '\Gamma_\rho';
            property_unit = '';
        case {'-gammas1'}
                property_name = 'Dimensionless slope RH';
                property_latex = '-\Gamma_\rho';
                property_unit = '';
        case {'gammas2', 'gammas'}
            property_name = 'Dimensionless slope RH';
            property_latex = '\Gamma';
            property_unit = '';
        case {'gammas3'}
            property_name = 'Dimensionless slope RH';
            property_latex = '\Gamma_p';
            property_unit = '';
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
