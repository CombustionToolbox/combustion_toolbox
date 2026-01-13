classdef ShockTurbulenceModelVorticalEntropic < combustiontoolbox.shockturbulence.ShockTurbulenceModel
    % The :mat:func:`ShockTurbulenceModelVorticalEntropic` class characterizes 
    % turbulence amplification across a shock wave interacting with weak turbulence
    % comprised of vortical-entropic disturbances using linear theory.
    %
    % These models are based on our previous theoretical works [1-2]
    % and have been extended to multi-component mixtures [3-6] using the
    % Combustion Toolbox [7-8].
    %
    % References:
    %     [1] Huete, C., Cuadra, A., Vera, M., Urzay, & J. (2021). Thermochemical
    %         effects on hypersonic shock waves interacting with weak turbulence.
    %         Physics of Fluids 33, 086111 (featured article). DOI: 10.1063/5.0059948.
    %
    %     [2] Huete, C., Velikovich, A. L., & Wouchuk, J. G. (2011). Analytical linear theory
    %         for the interaction of a planar shock wave with a two-or three-dimensional
    %         random isotropic density field. Physical Review E—Statistical, Nonlinear, and
    %         Soft Matter Physics, 83(5), 056320. DOI: 10.1103/PhysRevE.83.056320.
    %
    %     [3] Cuadra, A., Vera, M., Di Renzo, M., & Huete, C. (2023). Linear Theory
    %         of Hypersonic Shocks Interacting with Turbulence in Air. In 2023 AIAA
    %         SciTech Forum, National Harbor, USA. DOI: 10.2514/6.2023-0075.
    %
    %     [4] Cuadra, A., Williams, C. T., Di Renzo, M. & Huete, C. (2024). Compressibility
    %         and vibrational-excitation effects in hypersonic shock-turbulence interaction.
    %         Tech. Rep. Summer Program Proceedings, Center for Turbulence Research,
    %         Stanford University.
    %
    %     [5] Cuadra, A., Williams, C. T., Di Renzo, M., & Huete, C. The role of compressibility
    %         and vibrational-excitation in hypersonic shock–turbulence interactions.
    %         Journal of Fluid Mechanics (under review).
    %
    %     [6] Cuadra, A., Di Renzo, M., Hoste, J. J. O., Williams, C. T., Vera, M., & Huete, C. (2025).
    %         Review of shock-turbulence interaction with a focus on hypersonic flow. Physics of Fluids, 37(4).
    %         DOI: 10.1063/5.0255816.
    %
    %     [7] Cuadra, A., Huete, C., & Vera, M. (2026). Combustion Toolbox: An open-source
    %         thermochemical code for gas-and condensed-phase problems involving chemical equilibrium. 
    %         Computer Physics Communications 320, 110004. DOI:10.1016/j.cpc.2025.110004.
    %
    %     [8] Cuadra, A., Huete, C., Vera, M. (2022). Combustion Toolbox:
    %         A MATLAB-GUI based open-source tool for solving gaseous
    %         combustion problems. Zenodo. DOI: 10.5281/zenodo.5554911.
    
    properties
        typeChiFunction
    end

    properties (Access = private)
        chiFunction
    end

    methods

        function obj = ShockTurbulenceModelVorticalEntropic(varargin)
            % Constructor
            %
            % Optional Args:
            %   varargin (optional): key-value pairs to initialize the database
            %
            % Returns:
            %     obj (ShockTurbulenceModelVorticalEntropic): ShockTurbulenceModelVorticalEntropic object
            %
            % Examples:
            %     * shockTurbulenceModel = combustiontoolbox.shockturbulence.ShockTurbulenceModelVorticalEntropic();
            %     * shockTurbulenceModel = combustiontoolbox.shockturbulence.ShockTurbulenceModelVorticalEntropic('dimensions', 3);
            %     * shockTurbulenceModel = combustiontoolbox.shockturbulence.ShockTurbulenceModelVorticalEntropic('tolIntegralRelative', 1e-6);
            %     * shockTurbulenceModel = combustiontoolbox.shockturbulence.ShockTurbulenceModelVorticalEntropic('tolIntegralAbsolute', 1e-10);

            % Default
            defaultTypeChiFunction = 'chi';

            % Call superclass constructor
            obj@combustiontoolbox.shockturbulence.ShockTurbulenceModel('problemType', 'vortical_entropic', varargin{:});

            % Parse inputs
            p = inputParser;
            addParameter(p, 'typeChiFunction', defaultTypeChiFunction, @(x) ischar(x) && any(strcmpi(x, {'chi', 'scaleM1'})));
            parse(p, varargin{:});

            % Set properties
            obj.typeChiFunction = p.Results.typeChiFunction;

            % Set probability density function (PDF) for the model
            setPDF(obj);

            % Set chi function
            setChiFunction(obj);
        end

        function obj = setPDF(obj)
            % Define the probability density function (PDF) used to weigh different
            % wavenumber contributions in the post-shock turbulence model
            %
            % Args:
            %     obj (ShockTurbulenceModelVorticalEntropic): ShockTurbulenceModelVorticalEntropic object
            %
            % Returns:
            %     obj (ShockTurbulenceModelVorticalEntropic): Updated ShockTurbulenceModelVorticalEntropic object with PDF set
            
            switch obj.dimensions
                case 2
                    % 2D PDF
                    error('To be implemented.');
                case 3
                    % 3D PDF
                    obj.pdf = @(R, M2, zeta) 3/2 * ( (M2.^4 .* R.^4 .* sqrt(1 - M2.^2)) ./ (M2.^2 .* R.^2 + zeta.^2 .* (1 - M2.^2)).^(5/2) );
                otherwise
                    error('Invalid dimensions. Only 2D and 3D are supported.');
            end

        end

        function obj = setChiFunction(obj)
            % Set the chi function for the model.
            %
            % Args:
            %     obj (ShockTurbulenceModelVorticalEntropic): ShockTurbulenceModelVorticalEntropic object
            %
            % Returns:
            %     obj (ShockTurbulenceModelVorticalEntropic): Updated object with chi function set
            
            switch lower(obj.typeChiFunction)
                case 'chi'
                    % Chi function
                    obj.chiFunction = @(R, M2, beta, chi) chi;
                case 'scalem1'
                    % Scale M1 function
                    obj.chiFunction = @(R, M2, beta, chi) chi ./ combustiontoolbox.shockturbulence.ShockTurbulenceModelVorticalEntropic.M1(R, M2, beta);
            end
            
        end
            
        function [averages, mixArray1, mixArray2] = getAverages(obj, jumpConditions, mixArray1, mixArray2)
            % Compute the post-shock turbulence statistics for vortical-entropic disturbances
            %
            % Args:
            %     obj (ShockTurbulenceModelVorticalEntropic): ShockTurbulenceModelVorticalEntropic object
            %     jumpConditions (struct): Structure with jump conditions across the shock wave
            %     mixArray1 (Mixture): Pre-shock Mixture objects
            %     mixArray2 (Mixture): Post-shock Mixture objects
            %
            % Returns:
            %     Tuple containing:
            %
            %     * averages (struct): Structure with averages of post-shock turbulence statistics (e.g., Reynolds stresses, turbulent kinetic energy, enstrophy, etc.)
            %     * mixArray1 (Mixture): Pre-shock Mixture objects
            %     * mixArray2 (Mixture): Post-shock Mixture objects
            %
            % Example:
            %     [averages, mixArray1, mixArray2] = getAverages(ShockTurbulenceModelVorticalEntropic(), jumpConditions, mixArray1, mixArray2);

            % Set properties
            R = jumpConditions.Rratio;        % Density ratio (rho2/rho1)
            M2 = jumpConditions.M2;           % Post-shock Mach number
            Gammas = jumpConditions.Gammas2;  % Dimensionless slope of the Hugoniot curve (partial derivative at constant rho1, p1)
            Gammas1 = jumpConditions.Gammas1; % Dimensionless slope of the Hugoniot curve (partial derivative at constant rho2, p1)
            beta = jumpConditions.beta;       % Ratio of speed of sound across the shock
            chi = [mixArray1.chi];            % Mean correlation coefficient between the vortical and entropic modes

            % Definitions
            numCases = length(R);
            
            % Evaluate chi in acordance with the selected chi function
            chiEval = obj.chiFunction(R, M2, beta, chi);

            % Compute acoustic and vortical modes of the longitudinal and
            % transverse components of the turbulent kinetic energy (TKE)
            % amplification
            for i = numCases:-1:1
                averages.R11r(i) = R11r(obj, R(i), M2(i), Gammas(i), Gammas1(i), beta(i), chiEval(i));
                averages.R11a(i) = R11a(obj, R(i), M2(i), Gammas(i), Gammas1(i), beta(i), chiEval(i));
                averages.RTTr(i) = RTTr(obj, R(i), M2(i), Gammas(i), Gammas1(i), beta(i), chiEval(i));
                averages.RTTa(i) = RTTa(obj, R(i), M2(i), Gammas(i), Gammas1(i), beta(i), chiEval(i));
                averages.enstrophy33(i) = enstrophy33(obj, R(i), M2(i), Gammas(i), Gammas1(i), beta(i), chiEval(i));
            end

            % Compute longitudinal contribution of the turbulent kinetic energy (TKE) amplification
            averages.R11 = averages.R11r + averages.R11a;

            % Compute transverse contribution of the turbulent kinetic energy (TKE) amplification
            averages.RTT = averages.RTTr + averages.RTTa;

            % Compute amplification of the turbulent kinetic energy (TKE)
            averages.K = 1/3 * (averages.R11 + 2 * averages.RTT);
            averages.Ka = 1/3 * (averages.R11a + 2 * averages.RTTa);
            averages.Kr = 1/3 * (averages.R11r + 2 * averages.RTTr);
            
            % Compute enstrophy
            averages.enstrophyTT = 1/4 * (R.^2 + 3 * averages.enstrophy33);
            averages.enstrophy = 1/3 * (1 + 2 * averages.enstrophyTT);

            % Compute anisotropy
            averages.anisotropy = 1 - (4 * averages.R11) ./ (3 * averages.K + averages.R11);

            % Get Kolmogorov length scale ratio across the shock
            averages.kolmogorovLengthRatio = obj.getKolmogorovLength(averages, mixArray1, mixArray2);
            
            % Set the post-shock turbulence statistics in the Mixture array
            mixArray2 = obj.setAverages2MixArray(averages, mixArray2);
        end

    end

    methods (Access = private)

        function value = R11rl(obj, R, M2, Gammas, Gammas1, beta, chi)
            % Compute the rotational contribution of the longitudinal TKE amplification ratio (longwave: zeta < 1)
            fun = @(zeta) (obj.Delta_u_l1(R, M2, Gammas, Gammas1, beta, chi, zeta).^2 + obj.Delta_u_l2(R, M2, Gammas, Gammas1, beta, chi, zeta).^2) .* obj.pdf(R, M2, zeta); 
            value = beta.^2 .* obj.integrate(fun, 0, 1);
        end
        
        function value = R11rs(obj, R, M2, Gammas, Gammas1, beta, chi)
            % Compute the rotational contribution of the longitudinal TKE amplification ratio (shortwave: zeta > 1)
            fun = @(zeta) obj.Delta_u_s(R, M2, Gammas, Gammas1, beta, chi, zeta).^2 .* obj.pdf(R, M2, zeta);
            value = beta.^2 .* obj.integrate(fun, 1, Inf);
        end
        
        function value = R11r(obj, R, M2, Gammas, Gammas1, beta, chi)
            % Compute the total rotational contribution of the longitudinal TKE amplification ratio (longwave + shortwave)
            value = R11rl(obj, R, M2, Gammas, Gammas1, beta, chi) + R11rs(obj, R, M2, Gammas, Gammas1, beta, chi);
        end
        
        function value = R11a(obj, R, M2, Gammas, Gammas1, beta, chi)
            % Compute the acoustic contribution of the longitudinal TKE amplification ratio (shortwave: zeta > 1)
            fun = @(zeta) obj.Delta_ua(R, M2, Gammas, Gammas1, beta, chi, zeta).^2 .* obj.pdf(R, M2, zeta);
            value = beta.^2 .* obj.integrate(fun, 1, Inf);
        end
        
        function value = R11(obj, R, M2, Gammas, Gammas1, beta, chi)
            % Compute the total longitudinal TKE amplification ratio (rotational + acoustic)
            value = R11r(obj, R, M2, Gammas, Gammas1, beta, chi) + R11a(obj, R, M2, Gammas, Gammas1, beta, chi);
        end
        
        function value = RTTrl(obj, R, M2, Gammas, Gammas1, beta, chi)
            % Compute the rotational contribution of the transverse TKE amplification ratio (longwave: zeta < 1)
            fun = @(zeta) (obj.Delta_v_l1(R, M2, Gammas, Gammas1, beta, chi, zeta).^2 + obj.Delta_v_l2(R, M2, Gammas, Gammas1, beta, chi, zeta).^2 + 3/2 * beta^(-2)) .* obj.pdf(R, M2, zeta);
            value = 0.5 * beta.^2 .* obj.integrate(fun, 0, 1);
        end
        
        function value = RTTrs(obj, R, M2, Gammas, Gammas1, beta, chi)
            % Compute the rotational contribution of the transverse TKE amplification ratio (shortwave: zeta > 1)
            fun = @(zeta) (obj.Delta_v_s(R, M2, Gammas, Gammas1, beta, chi, zeta).^2 + 3/2 * beta^(-2)).* obj.pdf(R, M2, zeta);
            value = 0.5 * beta.^2 .* obj.integrate(fun, 1, Inf);
        end
        
        function value = RTTr(obj, R, M2, Gammas, Gammas1, beta, chi)
            % Compute the total rotational contribution of the transverse TKE amplification ratio (longwave + shortwave)
            value = RTTrl(obj, R, M2, Gammas, Gammas1, beta, chi) + RTTrs(obj, R, M2, Gammas, Gammas1, beta, chi);
        end
        
        function value = RTTa(obj, R, M2, Gammas, Gammas1, beta, chi)
            % Compute the acoustic contribution of the transverse TKE amplification ratio (shortwave: zeta > 1)
            fun = @(zeta) obj.Delta_va(R, M2, Gammas, Gammas1, beta, chi, zeta).^2 .* obj.pdf(R, M2, zeta);
            value =  0.5 * beta.^2 .* obj.integrate(fun, 1, Inf);
        end
        
        function value = RTT(obj, R, M2, Gammas, Gammas1, beta, chi)
            % Compute the total transverse TKE amplification ratio (rotational + acoustic)
            value = RTTr(obj, R, M2, Gammas, Gammas1, beta, chi) + RTTa(obj, R, M2, Gammas, Gammas1, beta, chi);
        end

        function value = enstrophy33l(obj, R, M2, Gammas, Gammas1, beta, chi)
            % Compute the longwave contribution of the enstrophy amplification ratio (longwave: zeta < 1)
            fun = @(zeta) (obj.Delta_Omega_l1(R, M2, Gammas, Gammas1, beta, chi, zeta).^2 + obj.Delta_Omega_l2(R, M2, Gammas, Gammas1, beta, chi, zeta).^2 ) .* sin(obj.thetaOfzeta(R, M2, zeta)).^2 .* obj.pdf(R, M2, zeta);
            value = 2/3 * beta.^2 .* obj.integrate(fun, 0, 1);
        end
        
        function value = enstrophy33s(obj, R, M2, Gammas, Gammas1, beta, chi)
            % Compute the shortwave contribution of the enstrophy amplification ratio (shortwave: zeta > 1)
            fun = @(zeta) obj.Delta_Omega_s(R, M2, Gammas, Gammas1, beta, chi, zeta).^2 .* sin(obj.thetaOfzeta(R, M2, zeta)).^2 .* obj.pdf(R, M2, zeta);
            value = 2/3 * beta.^2 .* obj.integrate(fun, 1, Inf);
        end

        function value = enstrophy33(obj, R, M2, Gammas, Gammas1, beta, chi)
            % Compute the total enstrophy amplification ratio (longwave + shortwave)
            value = enstrophy33l(obj, R, M2, Gammas, Gammas1, beta, chi) + enstrophy33s(obj, R, M2, Gammas, Gammas1, beta, chi);
        end
        
        function printCheck(obj, R, M2, Gammas, Gammas1, beta, chi, zeta)
            % Debugging function to print the results of the calculations
            fprintf('\nCHECK SOLUTION\n\n');
    
            % Printing sigma values
            fprintf('sigma_a        = %.6f\n', obj.sigma_a(R, M2, Gammas));
            fprintf('sigma_b        = %.6f\n', obj.sigma_b(M2, Gammas));
            fprintf('sigma_c        = %.6f\n', obj.sigma_c(R, M2, Gammas));
            fprintf('\n');
            
            % Printing pi values
            fprintf('pi_l1          = %.6f\n', obj.pi_l1(R, M2, Gammas, Gammas1, beta, chi, zeta));
            fprintf('pi_l2          = %.6f\n', obj.pi_l2(R, M2, Gammas, Gammas1, beta, chi, zeta));
            fprintf('pi_s           = %.6f\n', obj.pi_s(R, M2, Gammas, Gammas1, beta, chi, zeta));
            fprintf('ka             = %.6f\n', obj.ka(M2, zeta));
            fprintf('wa             = %.6f\n', obj.wa(M2, zeta));
            fprintf('\n');
            
            % Printing Delta values
            fprintf('Delta_ua       = %.6f\n', obj.Delta_ua(R, M2, Gammas, Gammas1, beta, chi, zeta));
            fprintf('Delta_va       = %.6f\n', obj.Delta_va(R, M2, Gammas, Gammas1, beta, chi, zeta));
            fprintf('\n');
            
            % Printing Omega values
            fprintf('Omega_1        = %.6f\n', obj.Omega_1(R, M2, beta, zeta));
            fprintf('Omega_2        = %.6f\n', obj.Omega_2(R, M2, Gammas));
            fprintf('Delta_Omega_l1 = %.6f\n', obj.Delta_Omega_l1(R, M2, Gammas, Gammas1, beta, chi, zeta));
            fprintf('Delta_Omega_l2 = %.6f\n', obj.Delta_Omega_l2(R, M2, Gammas, Gammas1, beta, chi, zeta));
            fprintf('Delta_Omega_s  = %.6f\n', obj.Delta_Omega_s(R, M2, Gammas, Gammas1, beta, chi, zeta));
            fprintf('\n');
            
            % General Delta
            fprintf('Delta          = %.6f\n', obj.Delta(M2, zeta));
            fprintf('\n');   
            
            % Printing Delta_u values
            fprintf('Delta_u_l1     = %.6f\n', obj.Delta_u_l1(R, M2, Gammas, Gammas1, beta, chi, zeta));
            fprintf('Delta_u_l2     = %.6f\n', obj.Delta_u_l2(R, M2, Gammas, Gammas1, beta, chi, zeta));
            fprintf('Delta_u_s      = %.6f\n', obj.Delta_u_s(R, M2, Gammas, Gammas1, beta, chi, zeta));
            fprintf('Delta_u        = %.6f\n', obj.Delta_u(R, M2, Gammas, Gammas1, beta, chi, zeta));
            fprintf('\n');
            
            % Printing Delta_v values
            fprintf('Delta_v_l1     = %.6f\n', obj.Delta_v_l1(R, M2, Gammas, Gammas1, beta, chi, zeta));
            fprintf('Delta_v_l2     = %.6f\n', obj.Delta_v_l2(R, M2, Gammas, Gammas1, beta, chi, zeta));
            fprintf('Delta_v_s      = %.6f\n', obj.Delta_v_s(R, M2, Gammas, Gammas1, beta, chi, zeta));
            fprintf('Delta_v        = %.6f\n', obj.Delta_v(R, M2, Gammas, Gammas1, beta, chi, zeta));
            fprintf('\n');
            
            % Longitudinal TKE amplification
            fprintf('R11rl          = %.6f\n', obj.R11rl(R, M2, Gammas, Gammas1, beta, chi));
            fprintf('R11rs          = %.6f\n', obj.R11rs(R, M2, Gammas, Gammas1, beta, chi));
            fprintf('R11r           = %.6f\n', obj.R11r(R, M2, Gammas, Gammas1, beta, chi));
            fprintf('R11a           = %.6f\n', obj.R11a(R, M2, Gammas, Gammas1, beta, chi));
            fprintf('R11            = %.6f\n', obj.R11(R, M2, Gammas, Gammas1, beta, chi));
            fprintf('\n');
            
            % Transverse TKE amplification
            fprintf('RTTrl          = %.6f\n', obj.RTTrl(R, M2, Gammas, Gammas1, beta, chi));
            fprintf('RTTrs          = %.6f\n', obj.RTTrs(R, M2, Gammas, Gammas1, beta, chi));
            fprintf('RTTr           = %.6f\n', obj.RTTr(R, M2, Gammas, Gammas1, beta, chi));
            fprintf('RTTa           = %.6f\n', obj.RTTa(R, M2, Gammas, Gammas1, beta, chi));
            fprintf('RTT            = %.6f\n', obj.RTT(R, M2, Gammas, Gammas1, beta, chi));
        end

    end

    methods (Static, Access = ?combustiontoolbox.shockturbulence)

        function value = M1(R, M2, beta)
            value = beta .* R .* M2;
        end
        
        function value = thetaCritical(R, M2)
            value = atan(M2 .* R ./ sqrt(1 - M2.^2)); 
        end
        
        function value = thetaOfzeta(R, M2, zeta)
            value = atan((M2 .* R) ./ (sqrt(1 - M2.^2) .* zeta));
        end
        
        function value = sigma_a(R, M2, Gammas)
            value =  (R  ./ (R  - 1)) .* ((1 - Gammas) ./ (2*M2));
        end
        
        function value = sigma_b(M2, Gammas)
            value =  (1 + Gammas) ./ (2*M2);
        end
        
        function value = sigma_c(R, M2, Gammas)
            % Definitions
            sigma_a = combustiontoolbox.shockturbulence.ShockTurbulenceModelVorticalEntropic.sigma_a(R, M2, Gammas);
    
            % Compute sigma_c
            value =  (((M2.^2 .* (R - 1))) ./ (1 - M2.^2)) .* sigma_a;
        end

        function value = sigma_d(R, M2, Gammas, Gammas1)
            % Compute sigma_d
            value = M2 .* R .* (Gammas + Gammas1) ./ (2 * Gammas1);
        end

        function value = alpha1(R, M2, beta, zeta)
            % Compute alpha1
            value = (R - 1) ./ (beta .* R) .* (zeta.^2 - (R .* M2.^2) ./ (1 - M2.^2));
        end

        function value = alpha2(R, M2, Gammas, Gammas1, zeta)
            % Definitions
            sigma_d = combustiontoolbox.shockturbulence.ShockTurbulenceModelVorticalEntropic.sigma_d(R, M2, Gammas, Gammas1);

            % Compute alpha2
            value = - sigma_d .* (zeta.^2 - (R .* M2.^2) ./ (1 - M2.^2) .* ((M2 .* (R - 1)) ./ sigma_d - 1));
        end
        
        function value = pi_l1(R, M2, Gammas, Gammas1, beta, chi, zeta)
            % Definitions
            alpha1 = combustiontoolbox.shockturbulence.ShockTurbulenceModelVorticalEntropic.alpha1(R, M2, beta, zeta);
            alpha2 = combustiontoolbox.shockturbulence.ShockTurbulenceModelVorticalEntropic.alpha2(R, M2, Gammas, Gammas1, zeta);
            sigma_b = combustiontoolbox.shockturbulence.ShockTurbulenceModelVorticalEntropic.sigma_b(M2, Gammas);
            sigma_c = combustiontoolbox.shockturbulence.ShockTurbulenceModelVorticalEntropic.sigma_c(R, M2, Gammas);

            % Compute pi_l1
            value = -((alpha1 + alpha2 .* chi) .* (sigma_b .* zeta.^2 - sigma_c)) ./ (zeta.^2 .* (1 - zeta.^2) + (sigma_b .* zeta.^2 - sigma_c).^2);
        end
        
        function value = pi_l2(R, M2, Gammas, Gammas1, beta, chi, zeta)
            % Definitions
            alpha1 = combustiontoolbox.shockturbulence.ShockTurbulenceModelVorticalEntropic.alpha1(R, M2, beta, zeta);
            alpha2 = combustiontoolbox.shockturbulence.ShockTurbulenceModelVorticalEntropic.alpha2(R, M2, Gammas, Gammas1, zeta);
            sigma_b = combustiontoolbox.shockturbulence.ShockTurbulenceModelVorticalEntropic.sigma_b(M2, Gammas);
            sigma_c = combustiontoolbox.shockturbulence.ShockTurbulenceModelVorticalEntropic.sigma_c(R, M2, Gammas);

            % Compute pi_l2
            value = ((alpha1 + alpha2 .* chi) .* (zeta .* sqrt(1 - zeta.^2))) ./ (zeta.^2 .* (1 - zeta.^2) + (sigma_b .* zeta.^2 - sigma_c).^2);
        end
        
        function value = pi_s(R, M2, Gammas, Gammas1, beta, chi, zeta)
            % Definitions
            alpha1 = combustiontoolbox.shockturbulence.ShockTurbulenceModelVorticalEntropic.alpha1(R, M2, beta, zeta);
            alpha2 = combustiontoolbox.shockturbulence.ShockTurbulenceModelVorticalEntropic.alpha2(R, M2, Gammas, Gammas1, zeta);
            sigma_b = combustiontoolbox.shockturbulence.ShockTurbulenceModelVorticalEntropic.sigma_b(M2, Gammas);
            sigma_c = combustiontoolbox.shockturbulence.ShockTurbulenceModelVorticalEntropic.sigma_c(R, M2, Gammas);

            % Compute pi_s
            value = -((alpha1 + alpha2 .* chi)) ./ (zeta .* sqrt(zeta.^2 - 1) + sigma_b .* zeta.^2 - sigma_c);
        end
        
        function value = ka(M2, zeta)
            value = (zeta .* M2  - sqrt(zeta.^2 - 1)) ./ (sqrt(1 - M2.^2));
        end
        
        function value = wa(M2, zeta)
            value = (zeta - M2  .* sqrt(zeta.^2 - 1)) ./ (sqrt(1 - M2.^2));
        end
        
        function value = Delta_ua(R, M2, Gammas, Gammas1, beta, chi, zeta)
            % Definitions
            ka = combustiontoolbox.shockturbulence.ShockTurbulenceModelVorticalEntropic.ka(M2, zeta);
            wa = combustiontoolbox.shockturbulence.ShockTurbulenceModelVorticalEntropic.wa(M2, zeta);
            pi_s = combustiontoolbox.shockturbulence.ShockTurbulenceModelVorticalEntropic.pi_s(R, M2, Gammas, Gammas1, beta, chi, zeta);

            % Compute Delta_ua value
            value = (ka ./ wa) .* pi_s;
        end
        
        function value = Delta_va(R, M2, Gammas, Gammas1, beta, chi, zeta)
            % Definitions
            wa = combustiontoolbox.shockturbulence.ShockTurbulenceModelVorticalEntropic.wa(M2, zeta);
            pi_s = combustiontoolbox.shockturbulence.ShockTurbulenceModelVorticalEntropic.pi_s(R, M2, Gammas, Gammas1, beta, chi, zeta);

            % Compute Delta_va
            value = (1 ./ wa) .* pi_s;
        end
        
        function value = Omega_1(R, M2, beta, zeta)
            value = R .* beta.^(-1) .* (1 + zeta.^2 .* ((1 - M2.^2) ./ (R.^2 .* M2.^2)));
        end
        
        function value = Omega_2(R, M2, Gammas)
            value = ((R - 1) .* (1 - Gammas)) ./ (2 * M2);
        end

        function value = Omega_3(R, M2, Gammas, Gammas1)
            % Definitions
            sigma_d = combustiontoolbox.shockturbulence.ShockTurbulenceModelVorticalEntropic.sigma_d(R, M2, Gammas, Gammas1);

            % Compute Omega_3
            value = -(R - 1) .* (R .* M2 - sigma_d);
        end
        
        function value = Delta_Omega_l1(R, M2, Gammas, Gammas1, beta, chi, zeta)
            % Definitions
            Omega_1 = combustiontoolbox.shockturbulence.ShockTurbulenceModelVorticalEntropic.Omega_1(R, M2, beta, zeta);
            Omega_2 = combustiontoolbox.shockturbulence.ShockTurbulenceModelVorticalEntropic.Omega_2(R, M2, Gammas);
            Omega_3 = combustiontoolbox.shockturbulence.ShockTurbulenceModelVorticalEntropic.Omega_3(R, M2, Gammas, Gammas1);
            pi_l1 = combustiontoolbox.shockturbulence.ShockTurbulenceModelVorticalEntropic.pi_l1(R, M2, Gammas, Gammas1, beta, chi, zeta);

            % Compute Delta_Omega_l1
            value = Omega_1 + Omega_2 .* pi_l1 + Omega_3 .* chi;
        end
        
        function value = Delta_Omega_l2(R, M2, Gammas, Gammas1, beta, chi, zeta)
            % Definitions
            Omega_2 = combustiontoolbox.shockturbulence.ShockTurbulenceModelVorticalEntropic.Omega_2(R, M2, Gammas);
            pi_l2 = combustiontoolbox.shockturbulence.ShockTurbulenceModelVorticalEntropic.pi_l2(R, M2, Gammas, Gammas1, beta, chi, zeta);

            % Compute Delta_Omega_l2
            value = Omega_2 .* pi_l2;
        end
        
        function value = Delta_Omega_s(R, M2, Gammas, Gammas1, beta, chi, zeta)
            % Definitions
            Omega_1 = combustiontoolbox.shockturbulence.ShockTurbulenceModelVorticalEntropic.Omega_1(R, M2, beta, zeta);
            Omega_2 = combustiontoolbox.shockturbulence.ShockTurbulenceModelVorticalEntropic.Omega_2(R, M2, Gammas);
            Omega_3 = combustiontoolbox.shockturbulence.ShockTurbulenceModelVorticalEntropic.Omega_3(R, M2, Gammas, Gammas1);
            pi_s = combustiontoolbox.shockturbulence.ShockTurbulenceModelVorticalEntropic.pi_s(R, M2, Gammas, Gammas1, beta, chi, zeta);

            % Compute Delta_Omega_s
            value = Omega_1 + Omega_2 .* pi_s + Omega_3 .* chi;
        end

        function value = Delta_Omega_l(R, M2, Gammas, Gammas1, beta, chi, zeta)
            % Definitions
            Delta_Omega_l1 = combustiontoolbox.shockturbulence.ShockTurbulenceModelVorticalEntropic.Delta_Omega_l1(R, M2, Gammas, Gammas1, beta, chi, zeta);
            Delta_Omega_l2 = combustiontoolbox.shockturbulence.ShockTurbulenceModelVorticalEntropic.Delta_Omega_l2(R, M2, Gammas, zeta);

            % Compute Delta_Omega_l
            value = sqrt(Delta_Omega_l1.^2 + Delta_Omega_l2.^2);
        end
        
        function value = Delta(M2, zeta)
            value = 1 + ((1 - M2.^2) ./ M2.^2) .* zeta.^2;
        end
        
        function value = Delta_u_l1(R, M2, Gammas, Gammas1, beta, chi, zeta)
            % Definitions
            Delta_Omega_l1 = combustiontoolbox.shockturbulence.ShockTurbulenceModelVorticalEntropic.Delta_Omega_l1(R, M2, Gammas, Gammas1, beta, chi, zeta);
            Delta = combustiontoolbox.shockturbulence.ShockTurbulenceModelVorticalEntropic.Delta(M2, zeta);

            % Compute Delta_u_l1
            value = Delta_Omega_l1 ./ Delta;
        end
        
        function value = Delta_u_l2(R, M2, Gammas, Gammas1, beta, chi, zeta)
            % Definitions
            Delta_Omega_l2 = combustiontoolbox.shockturbulence.ShockTurbulenceModelVorticalEntropic.Delta_Omega_l2(R, M2, Gammas, Gammas1, beta, chi, zeta);
            Delta = combustiontoolbox.shockturbulence.ShockTurbulenceModelVorticalEntropic.Delta(M2, zeta);

            % Compute Delta_u_l2
            value = Delta_Omega_l2 ./ Delta;
        end
        
        function value = Delta_u_s(R, M2, Gammas, Gammas1, beta, chi, zeta)
            % Definitions
            Delta_Omega_s = combustiontoolbox.shockturbulence.ShockTurbulenceModelVorticalEntropic.Delta_Omega_s(R, M2, Gammas, Gammas1, beta, chi, zeta);
            Delta = combustiontoolbox.shockturbulence.ShockTurbulenceModelVorticalEntropic.Delta(M2, zeta);

            % Compute Delta_u_s
            value = Delta_Omega_s ./ Delta;
        end
        
        function value = Delta_u(R, M2, Gammas, Gammas1, beta, chi, zeta)
            
            if zeta > 1
                % Compute Delta_u for shortwave
                value = combustiontoolbox.shockturbulence.ShockTurbulenceModelVorticalEntropic.Delta_u_s(R, M2, Gammas, Gammas1, beta, chi, zeta);
                return
            end
        
            % Definitions
            Delta_u_l1 = combustiontoolbox.shockturbulence.ShockTurbulenceModelVorticalEntropic.Delta_u_l1(R, M2, Gammas, Gammas1, beta, chi, zeta);
            Delta_u_l2 = combustiontoolbox.shockturbulence.ShockTurbulenceModelVorticalEntropic.Delta_u_l2(R, M2, Gammas, Gammas1, beta, chi, zeta);

            % Compute Delta_u for longwave
            value = sqrt(Delta_u_l1.^2 + Delta_u_l2.^2);
        end
        
        function value = Delta_v_l1(R, M2, Gammas, Gammas1, beta, chi, zeta)
            % Definitions
            Delta_Omega_l1 = combustiontoolbox.shockturbulence.ShockTurbulenceModelVorticalEntropic.Delta_Omega_l1(R, M2, Gammas, Gammas1, beta, chi, zeta);
            Delta = combustiontoolbox.shockturbulence.ShockTurbulenceModelVorticalEntropic.Delta(M2, zeta);

            % Compute Delta_v_l1
            value = zeta.* ((sqrt(1 - M2.^2)) ./ M2) .* (Delta_Omega_l1 ./ Delta);
        end
        
        function value = Delta_v_l2(R, M2, Gammas, Gammas1, beta, chi, zeta)
            % Definitions
            Delta_Omega_l2 = combustiontoolbox.shockturbulence.ShockTurbulenceModelVorticalEntropic.Delta_Omega_l2(R, M2, Gammas, Gammas1, beta, chi, zeta);
            Delta = combustiontoolbox.shockturbulence.ShockTurbulenceModelVorticalEntropic.Delta(M2, zeta);

            % Compute Delta_v_l2
            value = zeta.* ((sqrt(1 - M2.^2)) ./ M2) .* (Delta_Omega_l2 ./ Delta);
        end
        
        function value = Delta_v_s(R, M2, Gammas, Gammas1, beta, chi, zeta)
            % Definitions
            Delta_Omega_s = combustiontoolbox.shockturbulence.ShockTurbulenceModelVorticalEntropic.Delta_Omega_s(R, M2, Gammas, Gammas1, beta, chi, zeta);
            Delta = combustiontoolbox.shockturbulence.ShockTurbulenceModelVorticalEntropic.Delta(M2, zeta);

            % Compute Delta_v_s
            value = zeta.* ((sqrt(1 - M2.^2)) ./ M2) .* (Delta_Omega_s ./ Delta);
        end
        
        function value = Delta_v(R, M2, Gammas, Gammas1, beta, chi, zeta)
            
            if zeta > 1
                % Compute Delta_v for shortwave
                value = combustiontoolbox.shockturbulence.ShockTurbulenceModelVorticalEntropic.Delta_v_s(R, M2, Gammas, Gammas1, beta, chi, zeta);
                return
            end
        
            % Definitions
            Delta_v_l1 = combustiontoolbox.shockturbulence.ShockTurbulenceModelVorticalEntropic.Delta_v_l1(R, M2, Gammas, Gammas1, beta, chi, zeta);
            Delta_v_l2 = combustiontoolbox.shockturbulence.ShockTurbulenceModelVorticalEntropic.Delta_v_l2(R, M2, Gammas, Gammas1, beta, chi, zeta);

            % Compute Delta_v for longwave
            value = sqrt(Delta_v_l1.^2 + Delta_v_l2.^2);
        end

    end

end
