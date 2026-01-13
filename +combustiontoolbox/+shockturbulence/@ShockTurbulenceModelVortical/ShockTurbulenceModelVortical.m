classdef ShockTurbulenceModelVortical < combustiontoolbox.shockturbulence.ShockTurbulenceModel
    % The :mat:func:`ShockTurbulenceModelVortical` class characterizes 
    % turbulence amplification across a shock wave interacting with weak turbulence
    % comprised of vortical disturbances using linear theory.
    %
    % These models are based on our previous theoretical work [1]
    % and have been extended to multi-component mixtures [2-5] using the
    % Combustion Toolbox [6-7].
    %
    % References:
    %     [1] Huete, C., Cuadra, A., Vera, M., Urzay, & J. (2021). Thermochemical
    %         effects on hypersonic shock waves interacting with weak turbulence.
    %         Physics of Fluids 33, 086111 (featured article). DOI: 10.1063/5.0059948.
    %
    %     [2] Cuadra, A., Vera, M., Di Renzo, M., & Huete, C. (2023). Linear Theory
    %         of Hypersonic Shocks Interacting with Turbulence in Air. In 2023 AIAA
    %         SciTech Forum, National Harbor, USA. DOI: 10.2514/6.2023-0075.
    %
    %     [3] Cuadra, A., Williams, C. T., Di Renzo, M. & Huete, C. (2024). Compressibility
    %         and vibrational-excitation effects in hypersonic shock-turbulence interaction.
    %         Tech. Rep. Summer Program Proceedings, Center for Turbulence Research,
    %         Stanford University.
    %
    %     [4] Cuadra, A., Williams, C. T., Di Renzo, M., & Huete, C. The role of compressibility
    %         and vibrational-excitation in hypersonic shock–turbulence interactions.
    %         Journal of Fluid Mechanics (under review).
    %
    %     [5] Cuadra, A., Di Renzo, M., Hoste, J. J. O., Williams, C. T., Vera, M., & Huete, C. (2025).
    %         Review of shock-turbulence interaction with a focus on hypersonic flow. Physics of Fluids, 37(4).
    %         DOI: 10.1063/5.0255816.
    %
    %     [6] Cuadra, A., Huete, C., & Vera, M. (2026). Combustion Toolbox: An open-source
    %         thermochemical code for gas-and condensed-phase problems involving chemical equilibrium. 
    %         Computer Physics Communications 320, 110004. DOI:10.1016/j.cpc.2025.110004.
    %
    %     [7] Cuadra, A., Huete, C., Vera, M. (2022). Combustion Toolbox:
    %         A MATLAB-GUI based open-source tool for solving gaseous
    %         combustion problems. Zenodo. DOI: 10.5281/zenodo.5554911.
    
    methods

        function obj = ShockTurbulenceModelVortical(varargin)
            % Constructor
            %
            % Optional Args:
            %   varargin (optional): key-value pairs to initialize the database
            %
            % Returns:
            %     obj (ShockTurbulenceModelVortical): ShockTurbulenceModelVortical object
            %
            % Examples:
            %     * shockTurbulenceModel = combustiontoolbox.shockturbulence.ShockTurbulenceModelVortical();
            %     * shockTurbulenceModel = combustiontoolbox.shockturbulence.ShockTurbulenceModelVortical('dimensions', 3);
            %     * shockTurbulenceModel = combustiontoolbox.shockturbulence.ShockTurbulenceModelVortical('tolIntegralRelative', 1e-6);
            %     * shockTurbulenceModel = combustiontoolbox.shockturbulence.ShockTurbulenceModelVortical('tolIntegralAbsolute', 1e-10);

            % Call superclass constructor
            obj@combustiontoolbox.shockturbulence.ShockTurbulenceModel('problemType', 'vortical', varargin{:});

            % Set probability density function (PDF) for the model
            setPDF(obj);
        end

        function obj = setPDF(obj)
            % Define the probability density function (PDF) used to weigh different
            % wavenumber contributions in the post-shock turbulence model
            %
            % Args:
            %     obj (ShockTurbulenceModelVortical): ShockTurbulenceModelVortical object
            %
            % Returns:
            %     obj (ShockTurbulenceModelVortical): Updated ShockTurbulenceModelVortical object with PDF set
            
            switch obj.dimensions
                case 2
                    % 2D PDF
                    obj.pdf = @(R, M2, zeta) ( (M2.^3 .* R.^3 .* sqrt(1 - M2.^2)) ./ (M2.^2 .* R.^2 + zeta.^2 .* (1 - M2.^2)).^(2) );
                case 3
                    % 3D PDF
                    obj.pdf = @(R, M2, zeta) 3/2 * ( (M2.^4 .* R.^4 .* sqrt(1 - M2.^2)) ./ (M2.^2 .* R.^2 + zeta.^2 .* (1 - M2.^2)).^(5/2) );
                otherwise
                    error('Invalid dimensions. Only 2D and 3D are supported.');
            end

        end

        function [averages, mixArray1, mixArray2] = getAverages(obj, jumpConditions, mixArray1, mixArray2)
            % Compute the post-shock turbulence statistics for vortical disturbances
            %
            % Args:
            %     obj (ShockTurbulenceModelVortical): ShockTurbulenceModelVortical object
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
            %     [averages, mixArray1, mixArray2] = getAverages(ShockTurbulenceModelVortical(), jumpConditions, mixArray1, mixArray2);

            % Set properties
            R = jumpConditions.Rratio;       % Density ratio (rho2/rho1)
            M2 = jumpConditions.M2;          % Post-shock Mach number
            Gammas = jumpConditions.Gammas2; % Dimensionless slope of the Hugoniot curve (partial derivative at constant rho1, p1)

            % Definitions
            numCases = length(R);

            % Compute acoustic and vortical modes of the longitudinal and
            % transverse components of the turbulent kinetic energy (TKE)
            % amplification
            for i = numCases:-1:1
                averages.R11r(i) = R11r(obj, R(i), M2(i), Gammas(i));
                averages.R11a(i) = R11a(obj, R(i), M2(i), Gammas(i));
                averages.RTTr(i) = RTTr(obj, R(i), M2(i), Gammas(i));
                averages.RTTa(i) = RTTa(obj, R(i), M2(i), Gammas(i));
                averages.enstrophy33(i) = enstrophy33(obj, R(i), M2(i), Gammas(i));
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

        function value = R11rl(obj, R, M2, Gammas)
            % Compute the rotational contribution of the longitudinal TKE amplification ratio (longwave: zeta < 1)
            fun = @(zeta) (obj.Delta_u_l1(R, M2, Gammas, zeta).^2 + obj.Delta_u_l2(R, M2, Gammas, zeta).^2) .* obj.pdf(R, M2, zeta); 
            value = obj.integrate(fun, 0, 1);
        end
        
        function value = R11rs(obj, R, M2, Gammas)
            % Compute the rotational contribution of the longitudinal TKE amplification ratio (shortwave: zeta > 1)
            fun = @(zeta) obj.Delta_u_s(R, M2, Gammas, zeta).^2 .* obj.pdf(R, M2, zeta);
            value = obj.integrate(fun, 1, Inf);
        end
        
        function value = R11r(obj, R, M2, Gammas)
            % Compute the total rotational contribution of the longitudinal TKE amplification ratio (longwave + shortwave)
            value = R11rl(obj, R, M2, Gammas) + R11rs(obj, R, M2, Gammas);
        end
        
        function value = R11a(obj, R, M2, Gammas)
            % Compute the acoustic contribution of the longitudinal TKE amplification ratio (shortwave: zeta > 1)
            fun = @(zeta) obj.Delta_ua(R, M2, Gammas, zeta).^2 .* obj.pdf(R, M2, zeta);
            value = obj.integrate(fun, 1, Inf);
        end
        
        function value = R11(obj, R, M2, Gammas)
            % Compute the total longitudinal TKE amplification ratio (rotational + acoustic)
            value = R11r(obj, R, M2, Gammas) + R11a(obj, R, M2, Gammas);
        end
        
        function value = RTTrl(obj, R, M2, Gammas)
            % Compute the rotational contribution of the transverse TKE amplification ratio (longwave: zeta < 1)
            fun = @(zeta) (obj.Delta_v_l1(R, M2, Gammas, zeta).^2 + obj.Delta_v_l2(R, M2, Gammas, zeta).^2 + 3/2) .* obj.pdf(R, M2, zeta);
            value = 0.5 * obj.integrate(fun, 0, 1);
        end
        
        function value = RTTrs(obj, R, M2, Gammas)
            % Compute the rotational contribution of the transverse TKE amplification ratio (shortwave: zeta > 1)
            fun = @(zeta) (obj.Delta_v_s(R, M2, Gammas, zeta).^2 + 3/2).* obj.pdf(R, M2, zeta);
            value = 0.5 * obj.integrate(fun, 1, Inf);
        end
        
        function value = RTTr(obj, R, M2, Gammas)
            % Compute the total rotational contribution of the transverse TKE amplification ratio (longwave + shortwave)
            value = RTTrl(obj, R, M2, Gammas) + RTTrs(obj, R, M2, Gammas);
        end
        
        function value = RTTa(obj, R, M2, Gammas)
            % Compute the acoustic contribution of the transverse TKE amplification ratio (shortwave: zeta > 1)
            fun = @(zeta) obj.Delta_va(R, M2, Gammas, zeta).^2 .* obj.pdf(R, M2, zeta);
            value =  0.5 * obj.integrate(fun, 1, Inf);
        end
        
        function value = RTT(obj, R, M2, Gammas)
            % Compute the total transverse TKE amplification ratio (rotational + acoustic)
            value = RTTr(obj, R, M2, Gammas) + RTTa(obj, R, M2, Gammas);
        end

        function value = enstrophy33l(obj, R, M2, Gammas)
            % Compute the longwave contribution of the enstrophy amplification ratio (longwave: zeta < 1)
            fun = @(zeta) (obj.Delta_Omega_l1(R, M2, Gammas, zeta).^2 + obj.Delta_Omega_l2(R, M2, Gammas, zeta).^2 ) .* sin(obj.thetaOfzeta(R, M2, zeta)).^2 .* obj.pdf(R, M2, zeta);
            value = 2/3 * obj.integrate(fun, 0, 1);
        end
        
        function value = enstrophy33s(obj, R, M2, Gammas)
            % Compute the shortwave contribution of the enstrophy amplification ratio (shortwave: zeta > 1)
            fun = @(zeta) obj.Delta_Omega_s(R, M2, Gammas, zeta).^2 .* sin(obj.thetaOfzeta(R, M2, zeta)).^2 .* obj.pdf(R, M2, zeta);
            value = 2/3 * obj.integrate(fun, 1, Inf);
        end

        function value = enstrophy33(obj, R, M2, Gammas)
            % Compute the total enstrophy amplification ratio (longwave + shortwave)
            value = enstrophy33l(obj, R, M2, Gammas) + enstrophy33s(obj, R, M2, Gammas);
        end

        function printCheck(obj, R, M2, Gammas, zeta)
            % Debugging function to print the results of the calculations
            fprintf('\nCHECK SOLUTION\n\n');

            % Printing sigma values
            fprintf('sigma_a        = %.6f\n', obj.sigma_a(R, M2, Gammas));
            fprintf('sigma_b        = %.6f\n', obj.sigma_b(M2, Gammas));
            fprintf('sigma_c        = %.6f\n', obj.sigma_c(R, M2, Gammas));
            fprintf('\n')

            % Printing pi values
            fprintf('pi_l1          = %.6f\n', obj.pi_l1(R, M2, Gammas, zeta));
            fprintf('pi_l2          = %.6f\n', obj.pi_l2(R, M2, Gammas, zeta));
            fprintf('pi_s           = %.6f\n', obj.pi_s(R, M2, Gammas, zeta));
            fprintf('ka             = %.6f\n', obj.ka(M2, zeta));
            fprintf('wa             = %.6f\n', obj.wa(M2, zeta));
            fprintf('\n')

            % Printing Delta values
            fprintf('Delta_ua       = %.6f\n', obj.Delta_ua(R, M2, Gammas, zeta));
            fprintf('Delta_va       = %.6f\n', obj.Delta_va(R, M2, Gammas, zeta));
            fprintf('\n')

            % Printing Omega values
            fprintf('Omega_1        = %.6f\n', obj.Omega_1(R, M2, zeta));
            fprintf('Omega_2        = %.6f\n', obj.Omega_2(R, M2, Gammas));
            fprintf('Delta_Omega_l1 = %.6f\n', obj.Delta_Omega_l1(R, M2, Gammas, zeta));
            fprintf('Delta_Omega_l2 = %.6f\n', obj.Delta_Omega_l2(R, M2, Gammas, zeta));
            fprintf('Delta_Omega_s  = %.6f\n', obj.Delta_Omega_s(R, M2, Gammas, zeta));
            fprintf('\n')

            % General Delta
            fprintf('Delta          = %.6f\n', obj.Delta(M2, zeta));
            fprintf('\n')

            % Printing Delta_u values
            fprintf('Delta_u_l1     = %.6f\n', obj.Delta_u_l1(R, M2, Gammas, zeta));
            fprintf('Delta_u_l2     = %.6f\n', obj.Delta_u_l2(R, M2, Gammas, zeta));
            fprintf('Delta_u_s      = %.6f\n', obj.Delta_u_s(R, M2, Gammas, zeta));
            fprintf('Delta_u        = %.6f\n', obj.Delta_u(R, M2, Gammas, zeta));
            fprintf('\n')

            % Printing Delta_v values
            fprintf('Delta_v_l1     = %.6f\n', obj.Delta_v_l1(R, M2, Gammas, zeta));
            fprintf('Delta_v_l2     = %.6f\n', obj.Delta_v_l2(R, M2, Gammas, zeta));
            fprintf('Delta_v_s      = %.6f\n', obj.Delta_v_s(R, M2, Gammas, zeta));
            fprintf('Delta_v        = %.6f\n', obj.Delta_v(R, M2, Gammas, zeta));
            fprintf('\n')

            % Longitudinal TKE amplification
            fprintf('R11rl          = %.6f\n', obj.R11rl(R, M2, Gammas));
            fprintf('R11rs          = %.6f\n', obj.R11rs(R, M2, Gammas));
            fprintf('R11r           = %.6f\n', obj.R11r(R, M2, Gammas));
            fprintf('R11a           = %.6f\n', obj.R11a(R, M2, Gammas));
            fprintf('R11            = %.6f\n', obj.R11(R, M2, Gammas));
            fprintf('\n')

            % Transverse TKE amplification
            fprintf('RTTrl          = %.6f\n', obj.RTTrl(R, M2, Gammas));
            fprintf('RTTrs          = %.6f\n', obj.RTTrs(R, M2, Gammas));
            fprintf('RTTr           = %.6f\n', obj.RTTr(R, M2, Gammas));
            fprintf('RTTa           = %.6f\n', obj.RTTa(R, M2, Gammas));
            fprintf('RTT            = %.6f\n', obj.RTT(R, M2, Gammas));
        end

    end

    methods (Static, Access = ?combustiontoolbox.shockturbulence)

        function value = thetaCritical(R, M2)
            value = atan(M2 .* R ./ sqrt(1 - M2.^2)); 
        end
        
        function value = thetaOfzeta(R, M2, zeta)
            value = atan((M2 .* R) ./ (sqrt(1 - M2.^2) .* zeta));
        end

        function value = theta2zeta(R, M2, theta)
            value = M2 .* R ./ ( sqrt(1 - M2.^2) .* tan(theta) ) ; 
        end
        
        function value = wavenumber2zeta(R, M2, kParallel, kPerp)
            value = M2 .* R ./ sqrt(1 - M2.^2) .* kParallel ./ kPerp; 
        end
        
        function value = wavenumberRatioCritical(R, M2)
            value = 1 ./ ( M2 .* R ./ (sqrt(1 - M2.^2)) ); 
        end
        
        function value = kParallelCritical(R, M2, kPerp)
            value = kPerp ./ ( M2 .* R ./ (sqrt(1 - M2.^2)) ); 
        end
        
        function value = kPerpCritical(R, M2, kParallel)
            value = kParallel .* ( M2 .* R ./ (sqrt(1 - M2.^2)) ); 
        end
        
        function value = sigma_a(R, M2, Gammas)
            value =  (R  ./ (R  - 1)) .* ((1 - Gammas) ./ (2*M2));
        end
        
        function value = sigma_b(M2, Gammas)
            value =  (1 + Gammas) ./ (2*M2);
        end
        
        function value = sigma_c(R, M2, Gammas)
            % Definitions
            sigma_a = combustiontoolbox.shockturbulence.ShockTurbulenceModelVortical.sigma_a(R, M2, Gammas);
    
            % Compute sigma_c
            value =  (((M2.^2 .* (R - 1))) ./ (1 - M2.^2)) .* sigma_a;
        end
        
        function value = pi_l1(R, M2, Gammas, zeta)
            % Definitions
            sigma_b = combustiontoolbox.shockturbulence.ShockTurbulenceModelVortical.sigma_b(M2, Gammas);
            sigma_c = combustiontoolbox.shockturbulence.ShockTurbulenceModelVortical.sigma_c(R, M2, Gammas);

            % Compute pi_l1
            value = ((-(1 - R.^-1) .* (sigma_b .* zeta.^2 - sigma_c)) ./ (zeta.^2.* (1 - zeta.^2) + (sigma_b .* zeta.^2 - sigma_c).^2)) .* (zeta.^2 - (R .* M2.^2) ./ (1 - M2.^2));
        end
        
        function value = pi_l2(R, M2, Gammas, zeta)
            % Definitions
            sigma_b = combustiontoolbox.shockturbulence.ShockTurbulenceModelVortical.sigma_b(M2, Gammas);
            sigma_c = combustiontoolbox.shockturbulence.ShockTurbulenceModelVortical.sigma_c(R, M2, Gammas);

            % Compute pi_l2
            value = (((1 - R.^-1) .* zeta .* sqrt((1 - zeta.^2))) ./ (zeta.^2 .* (1 - zeta.^2) + (sigma_b .* zeta.^2 - sigma_c).^2)) .* (zeta.^2 - (R .* M2.^2) ./ (1 - M2.^2));
        end
        
        function value = pi_s(R, M2, Gammas, zeta)
            % Definitions
            sigma_b = combustiontoolbox.shockturbulence.ShockTurbulenceModelVortical.sigma_b(M2, Gammas);
            sigma_c = combustiontoolbox.shockturbulence.ShockTurbulenceModelVortical.sigma_c(R, M2, Gammas);

            % Compute pi_s
            value = ( (-(1 - R.^-1)) ./ (zeta .* sqrt(zeta.^2 - 1) + sigma_b .* zeta.^2 - sigma_c) ) .* ( zeta.^2 - (R .* M2.^2) ./ (1 - M2.^2) );
        end
        
        function value = ka(M2, zeta)
            value = (zeta .* M2  - sqrt(zeta.^2 - 1)) ./ (sqrt(1 - M2.^2));
        end
        
        function value = wa(M2, zeta)
            value = (zeta - M2  .* sqrt(zeta.^2 - 1)) ./ (sqrt(1 - M2.^2));
        end
        
        function value = Delta_ua(R, M2, Gammas, zeta)
            % Definitions
            ka = combustiontoolbox.shockturbulence.ShockTurbulenceModelVortical.ka(M2, zeta);
            wa = combustiontoolbox.shockturbulence.ShockTurbulenceModelVortical.wa(M2, zeta);
            pi_s = combustiontoolbox.shockturbulence.ShockTurbulenceModelVortical.pi_s(R, M2, Gammas, zeta);

            % Compute Delta_ua value
            value = (ka ./ wa) .* pi_s;
        end
        
        function value = Delta_va(R, M2, Gammas, zeta)
            % Definitions
            wa = combustiontoolbox.shockturbulence.ShockTurbulenceModelVortical.wa(M2, zeta);
            pi_s = combustiontoolbox.shockturbulence.ShockTurbulenceModelVortical.pi_s(R, M2, Gammas, zeta);

            % Compute Delta_va
            value = (1 ./ wa) .* pi_s;
        end
        
        function value = Omega_1(R, M2, zeta)
            value = R .* (1 + zeta.^2 .* ((1 - M2.^2) ./ (R.^2 .* M2.^2)));
        end
        
        function value = Omega_2(R, M2, Gammas)
            value = ((R - 1) .* (1 - Gammas)) ./ (2*M2);
        end
        
        function value = Delta_Omega_l1(R, M2, Gammas, zeta)
            % Definitions
            Omega_1 = combustiontoolbox.shockturbulence.ShockTurbulenceModelVortical.Omega_1(R, M2, zeta);
            Omega_2 = combustiontoolbox.shockturbulence.ShockTurbulenceModelVortical.Omega_2(R, M2, Gammas);
            pi_l1 = combustiontoolbox.shockturbulence.ShockTurbulenceModelVortical.pi_l1(R, M2, Gammas, zeta);

            % Compute Delta_Omega_l1
            value = Omega_2 .* pi_l1 + Omega_1;
        end
        
        function value = Delta_Omega_l2(R, M2, Gammas, zeta)
            % Definitions
            Omega_2 = combustiontoolbox.shockturbulence.ShockTurbulenceModelVortical.Omega_2(R, M2, Gammas);
            pi_l2 = combustiontoolbox.shockturbulence.ShockTurbulenceModelVortical.pi_l2(R, M2, Gammas, zeta);

            % Compute Delta_Omega_l2
            value = Omega_2 .* pi_l2;
        end
        
        function value = Delta_Omega_s(R, M2, Gammas, zeta)
            % Definitions
            Omega_1 = combustiontoolbox.shockturbulence.ShockTurbulenceModelVortical.Omega_1(R, M2, zeta);
            Omega_2 = combustiontoolbox.shockturbulence.ShockTurbulenceModelVortical.Omega_2(R, M2, Gammas);
            pi_s = combustiontoolbox.shockturbulence.ShockTurbulenceModelVortical.pi_s(R, M2, Gammas, zeta);

            % Compute Delta_Omega_s
            value = Omega_2 .* pi_s + Omega_1;
        end
        
        function value = Delta(M2, zeta)
            value = 1 + ((1 - M2.^2) ./ M2.^2) .* zeta.^2;
        end
        
        function value = Delta_u_l1(R, M2, Gammas, zeta)
            % Definitions
            Delta_Omega_l1 = combustiontoolbox.shockturbulence.ShockTurbulenceModelVortical.Delta_Omega_l1(R, M2, Gammas, zeta);
            Delta = combustiontoolbox.shockturbulence.ShockTurbulenceModelVortical.Delta(M2, zeta);

            % Compute Delta_u_l1
            value = Delta_Omega_l1 ./ Delta;
        end
        
        function value = Delta_u_l2(R, M2, Gammas, zeta)
            % Definitions
            Delta_Omega_l2 = combustiontoolbox.shockturbulence.ShockTurbulenceModelVortical.Delta_Omega_l2(R, M2, Gammas, zeta);
            Delta = combustiontoolbox.shockturbulence.ShockTurbulenceModelVortical.Delta(M2, zeta);

            % Compute Delta_u_l2
            value = Delta_Omega_l2 ./ Delta;
        end
        
        function value = Delta_u_s(R, M2, Gammas, zeta)
            % Definitions
            Delta_Omega_s = combustiontoolbox.shockturbulence.ShockTurbulenceModelVortical.Delta_Omega_s(R, M2, Gammas, zeta);
            Delta = combustiontoolbox.shockturbulence.ShockTurbulenceModelVortical.Delta(M2, zeta);

            % Compute Delta_u_s
            value = Delta_Omega_s ./ Delta;
        end
        
        function value = Delta_u(R, M2, Gammas, zeta)
            
            if zeta > 1
                % Compute Delta_u for shortwave
                value = combustiontoolbox.shockturbulence.ShockTurbulenceModelVortical.Delta_u_s(R, M2, Gammas, zeta);
                return
            end
        
            % Definitions
            Delta_u_l1 = combustiontoolbox.shockturbulence.ShockTurbulenceModelVortical.Delta_u_l1(R, M2, Gammas, zeta);
            Delta_u_l2 = combustiontoolbox.shockturbulence.ShockTurbulenceModelVortical.Delta_u_l2(R, M2, Gammas, zeta);

            % Compute Delta_u for longwave
            value = sqrt(Delta_u_l1.^2 + Delta_u_l2.^2);
        end
        
        function value = Delta_v_l1(R, M2, Gammas, zeta)
            % Definitions
            Delta_Omega_l1 = combustiontoolbox.shockturbulence.ShockTurbulenceModelVortical.Delta_Omega_l1(R, M2, Gammas, zeta);
            Delta = combustiontoolbox.shockturbulence.ShockTurbulenceModelVortical.Delta(M2, zeta);

            % Compute Delta_v_l1
            value = zeta.* ((sqrt(1 - M2.^2)) ./ M2) .* (Delta_Omega_l1 ./ Delta);
        end
        
        function value = Delta_v_l2(R, M2, Gammas, zeta)
            % Definitions
            Delta_Omega_l2 = combustiontoolbox.shockturbulence.ShockTurbulenceModelVortical.Delta_Omega_l2(R, M2, Gammas, zeta);
            Delta = combustiontoolbox.shockturbulence.ShockTurbulenceModelVortical.Delta(M2, zeta);

            % Compute Delta_v_l2
            value = zeta.* ((sqrt(1 - M2.^2)) ./ M2) .* (Delta_Omega_l2 ./ Delta);
        end
        
        function value = Delta_v_s(R, M2, Gammas, zeta)
            % Definitions
            Delta_Omega_s = combustiontoolbox.shockturbulence.ShockTurbulenceModelVortical.Delta_Omega_s(R, M2, Gammas, zeta);
            Delta = combustiontoolbox.shockturbulence.ShockTurbulenceModelVortical.Delta(M2, zeta);

            % Compute Delta_v_s
            value = zeta.* ((sqrt(1 - M2.^2)) ./ M2) .* (Delta_Omega_s ./ Delta);
        end
        
        function value = Delta_v(R, M2, Gammas, zeta)
            
            if zeta > 1
                % Compute Delta_v for shortwave
                value = combustiontoolbox.shockturbulence.ShockTurbulenceModelVortical.Delta_v_s(R, M2, Gammas, zeta);
                return
            end
        
            % Definitions
            Delta_v_l1 = combustiontoolbox.shockturbulence.ShockTurbulenceModelVortical.Delta_v_l1(R, M2, Gammas, zeta);
            Delta_v_l2 = combustiontoolbox.shockturbulence.ShockTurbulenceModelVortical.Delta_v_l2(R, M2, Gammas, zeta);

            % Compute Delta_v for longwave
            value = sqrt(Delta_v_l1.^2 + Delta_v_l2.^2);
        end

    end

end
