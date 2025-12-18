classdef ShockTurbulenceModelAcoustic < combustiontoolbox.shockturbulence.ShockTurbulenceModel
    % The :mat:class:`ShockTurbulenceModelAcoustic` class characterizes 
    % turbulence amplification across a shock wave interacting with weak turbulence
    % comprised of acoustic disturbances using linear theory.
    %
    % These models are based on our previous theoretical work [1]
    % and have been extended to multi-component mixtures [2] using the
    % Combustion Toolbox [3].
    %
    % References:
    %     [1] Huete, C., Wouchuk, J. G., & Velikovich, A. L. (2012). Analytical linear theory
    %         for the interaction of a planar shock wave with a two-or three-dimensional random
    %         isotropic acoustic wave field. Physical Review E—Statistical, Nonlinear, and Soft
    %         Matter Physics, 85(2), 026312. DOI: 10.1063/5.0059948.
    %
    %     [2] Cuadra, A., Williams, C. T., Di Renzo, M. & Huete, C. (2024). Compressibility
    %         and vibrational-excitation effects in hypersonic shock-turbulence interaction.
    %         Tech. Rep. Summer Program Proceedings, Center for Turbulence Research,
    %         Stanford University.
    %
    %     [3] Cuadra, A., Huete, C., Vera, M. (2022). Combustion Toolbox:
    %         A MATLAB-GUI based open-source tool for solving gaseous
    %         combustion problems. Zenodo. DOI: 10.5281/zenodo.5554911.
    
    methods

        function obj = ShockTurbulenceModelAcoustic(varargin)
            % Constructor
            %
            % Optional Args:
            %   varargin (optional): key-value pairs to initialize the database
            %
            % Returns:
            %     obj (ShockTurbulenceModelAcoustic): ShockTurbulenceModelAcoustic object
            %
            % Examples:
            %     * shockTurbulenceModel = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic();
            %     * shockTurbulenceModel = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic('dimensions', 3);
            %     * shockTurbulenceModel = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic('tolIntegralRelative', 1e-6);
            %     * shockTurbulenceModel = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic('tolIntegralAbsolute', 1e-10);

            % Call superclass constructor
            obj@combustiontoolbox.shockturbulence.ShockTurbulenceModel('problemType', 'acoustic', varargin{:});

            % Set probability density function (PDF) for the model
            setPDF(obj);
        end

        function obj = setPDF(obj)
            % Define the probability density function (PDF) used to weigh different
            % wavenumber contributions in the post-shock turbulence model
            %
            % Args:
            %     obj (ShockTurbulenceModelAcoustic): ShockTurbulenceModelAcoustic object
            %
            % Returns:
            %     obj (ShockTurbulenceModelAcoustic): Updated ShockTurbulenceModelAcoustic object with PDF set
            
            switch obj.dimensions
                case 2
                    % 2D PDF
                    obj.pdf = @(R, M2, beta, zeta) ( - 2 * ((M2 .* R .* beta).^2 - 1) .* M2 .* obj.sintheta(R, M2, beta, zeta).^2 ) ./ ( pi * (M2 .* R .* beta) .* sqrt(1 - M2.^2) .* ( obj.costheta(R, M2, beta, zeta) - M2 .* R .* beta) );
                case 3
                    % 3D PDF
                    obj.pdf = @(R, M2, beta, zeta) ( -((M2 .* R .* beta).^2 - 1) .* M2 .* obj.sintheta(R, M2, beta, zeta).^3 ) ./ ( (M2 .* R .* beta) .* sqrt(1 - M2.^2) .* ( obj.costheta(R, M2, beta, zeta) - M2 .* R .* beta) );
                otherwise
                    error('Invalid dimensions. Only 2D and 3D are supported.');
            end

        end
            
        function averages = getAverages(obj, R, M2, Gammas, Gammas1, Gammas3, beta)
            % getAverages Compute the amplification ratios.
            %
            %   results = getAverages(R, M2, Gammas, Gammas1, Gammas3, beta)
            %
            % Args:
            %     obj (ShockTurbulenceModelAcoustic): ShockTurbulenceModelAcoustic object
            %     R (float): Density ratio rho_2 / rho_1
            %     M2 (float): Post-shock Mach number
            %     Gammas (float): Inverse normalized Hugonito slope, see Eq. (22) in [1] (Eq. (4) in [2])
            %     Gammas3 (float): Inverse normalized Hugonito slope
            %     beta (float): Ratio of speed of sound in the post-shock state to the pre-shock state
            %
            % Returns:
            %     averages (struct): Structure with computed averages (K, L, T, etc.)
            %
            % Example:
            %     averages = getAverages(ShockTurbulenceModelAcoustic(), R, M2, Gammas, Gammas3, beta);
            
            % Parse input arguments
            p = inputParser;
            addRequired(p, 'R', @(x) isnumeric(x));
            addRequired(p, 'M2', @(x) isnumeric(x));
            addRequired(p, 'Gammas', @(x) isnumeric(x));
            addRequired(p, 'Gammas1', @(x) isnumeric(x));
            addRequired(p, 'Gammas3', @(x) isnumeric(x));
            addRequired(p, 'beta', @(x) isnumeric(x));
            parse(p, R, M2, Gammas, Gammas1, Gammas3, beta);

            % Set properties
            R = p.Results.R;
            M2 = p.Results.M2;
            Gammas = p.Results.Gammas;
            Gammas1 = p.Results.Gammas1;
            Gammas3 = p.Results.Gammas3;
            beta = p.Results.beta;

            % Definitions
            N = length(R);
            
            % Compute acoustic and vortical modes of the longitudinal and
            % transverse components of the turbulent kinetic energy (TKE)
            % amplification
            for i = N:-1:1
                averages.R11rl(i) = R11rl(obj, R(i), M2(i), Gammas(i), Gammas1(i), Gammas3(i), beta(i));
                averages.R11rs(i) = R11rs(obj, R(i), M2(i), Gammas(i), Gammas1(i), Gammas3(i), beta(i));
                averages.R11a(i)  = R11a(obj, R(i), M2(i), Gammas(i), Gammas1(i), Gammas3(i), beta(i));
                averages.RTTrl(i) = RTTrl(obj, R(i), M2(i), Gammas(i), Gammas1(i), Gammas3(i), beta(i));
                averages.RTTrs(i) = RTTrs(obj, R(i), M2(i), Gammas(i), Gammas1(i), Gammas3(i), beta(i));
                averages.RTTa(i)  = RTTa(obj, R(i), M2(i), Gammas(i), Gammas1(i), Gammas3(i), beta(i));
                averages.Krl(i)   = Krl(obj, R(i), M2(i), Gammas(i), Gammas1(i), Gammas3(i), beta(i));
                averages.Krs(i)   = Krs(obj, R(i), M2(i), Gammas(i), Gammas1(i), Gammas3(i), beta(i));
                averages.Ka(i)    = Ka(obj, R(i), M2(i), Gammas(i), Gammas1(i), Gammas3(i), beta(i));
                averages.enstrophyl(i) = enstrophyl(obj, R(i), M2(i), Gammas(i), Gammas1(i), Gammas3(i), beta(i));
                averages.enstrophys(i) = enstrophys(obj, R(i), M2(i), Gammas(i), Gammas1(i), Gammas3(i), beta(i));
            end

            % Compute longitudinal contribution of the turbulent kinetic energy (TKE) amplification
            averages.R11 = 3 * (averages.R11rl + averages.R11rs + averages.R11a);
            averages.R11r = 3 * (averages.R11rl + averages.R11rs);
            
            % Compute transverse contribution of the turbulent kinetic energy (TKE) amplification
            averages.RTT = 3/2 * (averages.RTTrl + averages.RTTrs + averages.RTTa);
            averages.RTTr = 3/2 * (averages.RTTrl + averages.RTTrs);

            % Compute amplification of the turbulent kinetic energy (TKE)
            averages.K = averages.Krl + averages.Krs + averages.Ka;
            averages.Kr = averages.Krl + averages.Krs;

            % Compute enstrophy amplification
            averages.enstrophy = averages.enstrophyl + averages.enstrophys;
            averages.enstrophyTT = 0; % There is not transverse contribution to enstrophy in acoustic case
            
            % Compute anisotropy
            averages.anisotropy = 1 - (4 * averages.R11) ./ (averages.K + averages.R11);

            % Return the averages structure as output.
            obj.averages = averages;
        end

    end

    methods (Access = private)

        function value = Krl(obj, R, M2, Gammas, Gammas1, Gammas3, beta)
            % Compute the rotational contribution of the longitudinal TKE amplification ratio (longwave: zeta < 1)
            fun = @(zeta) (1 + zeta.^2 ./ sinh( obj.thetas(M2) ).^2) .* obj.Qlrot(R, M2, Gammas, Gammas1, Gammas3, beta, zeta).^2 .* abs( obj.pdf(R, M2, beta, zeta) );
            value = 0.5 * beta.^2 * (obj.integrate(fun, -1, 0) + obj.integrate(fun, 0, 1));
        end
        
        function value = Krs(obj, R, M2, Gammas, Gammas1, Gammas3, beta)
            % Compute the rotational contribution of the longitudinal TKE amplification ratio (shortwave: zeta > 1)
            fun = @(zeta) (1 + zeta.^2 ./ sinh( obj.thetas(M2) ).^2) .* obj.Qsrot(R, M2, Gammas, Gammas1, Gammas3, beta, zeta).^2 .* abs( obj.pdf(R, M2, beta, zeta) );
            value = 0.5 * beta.^2 * (obj.integrate(fun, -Inf, -1) + obj.integrate(fun, 1, Inf));
        end
        
        function value = Ka(obj, R, M2, Gammas, Gammas1, Gammas3, beta)
            % Compute the acoustic contribution of the longitudinal TKE amplification ratio (shortwave: zeta > 1)
            fun_ne0 = @(zeta) obj.ne0(R, M2, Gammas, Gammas1, Gammas3, beta, zeta).^2 .* abs( obj.pdf(R, M2, beta, zeta) );
            fun_e0 = @(zeta) obj.e0(R, M2, Gammas, Gammas1, Gammas3, beta, zeta).^2 .* abs( obj.pdf(R, M2, beta, zeta) );
            value = 0.5 * beta.^2 * (obj.integrate(fun_ne0, -Inf, -1) + obj.integrate(fun_e0, 1, Inf));
        end
        
        function value = K(obj, R, M2, Gammas, Gammas1, Gammas3, beta)
            % Compute the total longitudinal TKE amplification ratio (rotational + acoustic)
            value = Krl(obj, R, M2, Gammas, Gammas1, Gammas3, beta) + Krs(obj, R, M2, Gammas, Gammas1, Gammas3, beta) + Ka(obj, R, M2, Gammas, Gammas1, Gammas3, beta);
        end
        
        function value = R11rl(obj, R, M2, Gammas, Gammas1, Gammas3, beta)
            % Compute the rotational contribution of the longitudinal TKE amplification ratio (longwave: zeta < 1)
            fun = @(zeta) obj.Qlrot(R, M2, Gammas, Gammas1, Gammas3, beta, zeta).^2 .* abs( obj.pdf(R, M2, beta, zeta) );
            value = 0.5 * beta.^2 * (obj.integrate(fun, -1, 0) + obj.integrate(fun, 0, 1));
        end
        
        function value = R11rs(obj, R, M2, Gammas, Gammas1, Gammas3, beta)
            % Compute the rotational contribution of the longitudinal TKE amplification ratio (shortwave: zeta > 1)
            fun = @(zeta) obj.Qsrot(R, M2, Gammas, Gammas1, Gammas3, beta, zeta).^2 .* abs( obj.pdf(R, M2, beta, zeta) );
            value = 0.5 * beta.^2 * (obj.integrate(fun, -Inf, -1) + obj.integrate(fun, 1, Inf));
        end
        
        function value = R11a(obj, R, M2, Gammas, Gammas1, Gammas3, beta)
            % Compute the acoustic contribution of the longitudinal TKE amplification ratio (shortwave: zeta > 1)
            fun = @(zeta) obj.Qax(R, M2, Gammas, Gammas1, Gammas3, beta, zeta).^2 .* abs( obj.pdf(R, M2, beta, zeta) );
            value = 0.5 * beta.^2 * (obj.integrate(fun, -Inf, -1) + obj.integrate(fun, 1, Inf));
        end
        
        function value = R11(obj, R, M2, Gammas, Gammas1, Gammas3, beta)
            % Compute the total longitudinal TKE amplification ratio (rotational + acoustic)
            value = 3 * (R11rl(obj, R, M2, Gammas, Gammas1, Gammas3, beta) + R11rs(obj, R, M2, Gammas, Gammas1, Gammas3, beta) + R11a(obj, R, M2, Gammas, Gammas1, Gammas3, beta));
        end
        
        function value = RTTrl(obj, R, M2, Gammas, Gammas1, Gammas3, beta)
            % Compute the rotational contribution of the transverse TKE amplification ratio (longwave: zeta < 1)
            fun = @(zeta) (zeta.^2 ./ sinh( obj.thetas(M2) ).^2) .* obj.Qlrot(R, M2, Gammas, Gammas1, Gammas3, beta, zeta).^2 .* abs( obj.pdf(R, M2, beta, zeta) );
            value = 0.5 * beta.^2 * (obj.integrate(fun, -1, 0) + obj.integrate(fun, 0, 1));
        end
        
        function value = RTTrs(obj, R, M2, Gammas, Gammas1, Gammas3, beta)
            % Compute the rotational contribution of the transverse TKE amplification ratio (shortwave: zeta > 1)
            fun = @(zeta) (zeta.^2 ./ sinh( obj.thetas(M2) ).^2) .* obj.Qsrot(R, M2, Gammas, Gammas1, Gammas3, beta, zeta).^2 .* abs( obj.pdf(R, M2, beta, zeta) );
            value = 0.5 * beta.^2 * (obj.integrate(fun, -Inf, -1) + obj.integrate(fun, 1, Inf));
        end
        
        function value = RTTa(obj, R, M2, Gammas, Gammas1, Gammas3, beta)
            % Compute the acoustic contribution of the transverse TKE amplification ratio (shortwave: zeta > 1)
            fun = @(zeta) (sqrt(1 - M2^2) ./ (M2 * abs(zeta) - sqrt(zeta.^2 - 1))).^2 .* obj.Qax(R, M2, Gammas, Gammas1, Gammas3, beta, zeta).^2 .* abs( obj.pdf(R, M2, beta, zeta) );
            value = 0.5 * beta.^2 * (obj.integrate(fun, -Inf, -1) + obj.integrate(fun, 1, Inf));
        end
        
        function value = RTT(obj, R, M2, Gammas, Gammas1, Gammas3, beta)
            % Compute the total transverse TKE amplification ratio (rotational + acoustic)
            value = 1.5 * (RTTrl(obj, R, M2, Gammas, Gammas1, Gammas3, beta) + RTTrs(obj, R, M2, Gammas, Gammas1, Gammas3, beta) + RTTa(obj, R, M2, Gammas, Gammas1, Gammas3, beta));
        end

        function value = enstrophyl(obj, R, M2, Gammas, Gammas1, Gammas3, beta)
            % Compute the longwave contribution of the enstrophy amplification ratio (longwave: zeta < 1)
            fun = @(zeta) ( (obj.b01(R, M2, Gammas, Gammas1, Gammas3, beta, zeta) - 1).^2 + obj.b02(R, M2, Gammas, Gammas1, Gammas3, beta, zeta).^2 ) .* obj.sintheta(R, M2, beta, zeta).^2 .* abs( obj.pdf(R, M2, beta, zeta) );
            value = 0.5 * obj.Omega_2(R, M2, Gammas).^2 .* (obj.integrate(fun, -1, 0) + obj.integrate(fun, 0, 1));
        end
        
        function value = enstrophys(obj, R, M2, Gammas, Gammas1, Gammas3, beta)
            % Compute the shortwave contribution of the enstrophy amplification ratio (shortwave: zeta > 1)
            funPos = @(zeta) ( obj.pe0(R, M2, Gammas, Gammas1, Gammas3, beta, zeta) - 1).^2 .* obj.sintheta(R, M2, beta, zeta).^2 .* abs( obj.pdf(R, M2, beta, zeta) );
            funNeg = @(zeta) ( obj.ne0(R, M2, Gammas, Gammas1, Gammas3, beta, zeta) - 1).^2 .* obj.sintheta(R, M2, beta, zeta).^2 .* abs( obj.pdf(R, M2, beta, zeta) );
            value = 0.5 * obj.Omega_2(R, M2, Gammas).^2 .* (obj.integrate(funNeg, -Inf, -1) + obj.integrate(funPos, 1, Inf));
        end

        function value = enstrophy(obj, R, M2, Gammas, Gammas1, Gammas3, beta)
            % Compute the total enstrophy amplification ratio (longwave + shortwave)
            value = enstrophyl(obj, R, M2, Gammas, Gammas1, Gammas3, beta) + enstrophys(obj, R, M2, Gammas, Gammas1, Gammas3, beta);
        end
        
        function printCheck(obj, R, M2, Gammas, Gammas1, Gammas3, beta, zeta)
            % Debugging function to print the results of the calculations
            fprintf('\nCHECK SOLUTION\n\n');

            % Sigma values
            fprintf('sigma_a        = %.6f\n', obj.sigma_a(R, M2, Gammas));
            fprintf('sigma_b        = %.6f\n', obj.sigma_b(M2, Gammas));
            fprintf('sigma_c        = %.6f\n', obj.sigma_c(R, M2, Gammas));
            fprintf('sigma_d        = %.6f\n', obj.sigma_d(R, M2, Gammas, Gammas1));
            fprintf('sigma_e        = %.6f\n', obj.sigma_e(R, Gammas, Gammas3, beta));
            fprintf('\n');
            
            % Pi values
            fprintf('pi_l1          = %.6f\n', obj.pi_l1(R, M2, Gammas, Gammas1, Gammas3, beta, zeta));
            fprintf('pi_l2          = %.6f\n', obj.pi_l2(R, M2, Gammas, Gammas1, Gammas3, beta, zeta));
            fprintf('pi_s           = %.6f\n', obj.pi_s(R, M2, Gammas, Gammas1, Gammas3, beta, zeta));
            fprintf('ka             = %.6f\n', obj.ka(M2, zeta));
            fprintf('wa             = %.6f\n', obj.wa(M2, zeta));
            fprintf('\n');
            
            % Delta_u and Delta_v values
            fprintf('Delta_ua       = %.6f\n', obj.Delta_ua(R, M2, Gammas, Gammas1, Gammas3, beta, zeta));
            fprintf('Delta_va       = %.6f\n', obj.Delta_va(R, M2, Gammas, Gammas1, Gammas3, beta, zeta));
            fprintf('\n');
            
            % Omega and Delta_Omega values
            fprintf('Omega_1        = %.6f\n', obj.Omega_1(R, M2, beta, zeta));
            fprintf('Omega_2        = %.6f\n', obj.Omega_2(R, M2, Gammas));
            fprintf('Delta_Omega_l1 = %.6f\n', obj.Delta_Omega_l1(R, M2, Gammas, Gammas1, Gammas3, beta, zeta));
            fprintf('Delta_Omega_l2 = %.6f\n', obj.Delta_Omega_l2(R, M2, Gammas, Gammas1, Gammas3, beta, zeta));
            fprintf('Delta_Omega_s  = %.6f\n', obj.Delta_Omega_s(R, M2, Gammas, Gammas1, Gammas3, beta, zeta));
            fprintf('\n');
            
            % Delta values
            fprintf('Delta          = %.6f\n', obj.Delta(M2, zeta));
            fprintf('\n');
            
            % Delta_u values
            fprintf('Delta_u_l1     = %.6f\n', obj.Delta_u_l1(R, M2, Gammas, Gammas1, Gammas3, beta, zeta));
            fprintf('Delta_u_l2     = %.6f\n', obj.Delta_u_l2(R, M2, Gammas, Gammas1, Gammas3, beta, zeta));
            fprintf('Delta_u_s      = %.6f\n', obj.Delta_u_s(R, M2, Gammas, Gammas1, Gammas3, beta, zeta));
            fprintf('Delta_u        = %.6f\n', obj.Delta_u(R, M2, Gammas, Gammas1, Gammas3, beta, zeta));
            fprintf('\n');
            
            % Delta_v values
            fprintf('Delta_v_l1     = %.6f\n', obj.Delta_v_l1(R, M2, Gammas, Gammas1, Gammas3, beta, zeta));
            fprintf('Delta_v_l2     = %.6f\n', obj.Delta_v_l2(R, M2, Gammas, Gammas1, Gammas3, beta, zeta));
            fprintf('Delta_v_s      = %.6f\n', obj.Delta_v_s(R, M2, Gammas, Gammas1, Gammas3, beta, zeta));
            fprintf('Delta_v        = %.6f\n', obj.Delta_v(R, M2, Gammas, Gammas1, Gammas3, beta, zeta));
            fprintf('\n');

            % e0 and ne0 values
            fprintf('e0             = %.6f\n', obj.e0(R, M2, Gammas, Gammas1, Gammas3, beta, zeta));
            fprintf('ne0            = %.6f\n', obj.ne0(R, M2, Gammas, Gammas1, Gammas3, beta, zeta));
            fprintf('\n');

            % Qsrot, Qlrot and Qax values
            fprintf('Qsrot          = %.6f\n', obj.Qsrot(R, M2, Gammas, Gammas1, Gammas3, beta, zeta));
            fprintf('Qlrot          = %.6f\n', obj.Qlrot(R, M2, Gammas, Gammas1, Gammas3, beta, zeta));
            fprintf('Qax            = %.6f\n', obj.Qax(R, M2, Gammas, Gammas1, Gammas3, beta, zeta));
            fprintf('\n');

            % Kl, Ks and Ka values
            fprintf('Kl             = %.6f\n', obj.Kl(R, M2, Gammas, Gammas1, Gammas3, beta));
            fprintf('Ks             = %.6f\n', obj.Ks(R, M2, Gammas, Gammas1, Gammas3, beta));
            fprintf('Ka             = %.6f\n', obj.Ka(R, M2, Gammas, Gammas1, Gammas3, beta));
            fprintf('\n');

            % K values
            fprintf('K              = %.6f\n', obj.K(R, M2, Gammas, Gammas1, Gammas3, beta));
            fprintf('\n');

            % R11l, R11s and R11a values
            fprintf('R11l           = %.6f\n', obj.R11l(R, M2, Gammas, Gammas1, Gammas3, beta));
            fprintf('R11s           = %.6f\n', obj.R11s(R, M2, Gammas, Gammas1, Gammas3, beta));
            fprintf('R11a           = %.6f\n', obj.R11a(R, M2, Gammas, Gammas1, Gammas3, beta));
            fprintf('\n');

            % R11 values
            fprintf('R11            = %.6f\n', obj.R11(R, M2, Gammas, Gammas1, Gammas3, beta));
            fprintf('\n');

            % RTTl, RTTs and RTTa values
            fprintf('RTTl           = %.6f\n', obj.RTTl(R, M2, Gammas, Gammas1, Gammas3, beta));
            fprintf('RTTs           = %.6f\n', obj.RTTs(R, M2, Gammas, Gammas1, Gammas3, beta));
            fprintf('RTTa           = %.6f\n', obj.RTTa(R, M2, Gammas, Gammas1, Gammas3, beta));
            fprintf('\n');

            % RTT values
            fprintf('RTT            = %.6f\n', obj.RTT(R, M2, Gammas, Gammas1, Gammas3, beta));
            fprintf('\n');
        end

    end

    methods (Static, Access = ?combustiontoolbox.shockturbulence)

        function value = M1(R, M2, beta)
            value = beta .* R .* M2;
        end
        
        function value = thetas(M2)
            value = atanh(M2);
        end
        
        function value = costheta(R, M2, beta, zeta)
            value = R .* M2 ./ (beta .* (R.^2 .* M2.^2 + (1 - M2.^2) .* zeta.^2)) .* (1 + zeta .* beta .* sqrt(1 - M2.^2) .* sqrt( (1 - M2.^2) ./  (R.^2 .* M2.^2) .* zeta.^2 + ((R .* beta .* M2).^2 - 1) ./ (R .* beta .* M2).^2));
        end
        
        function value = sintheta(R, M2, beta, zeta)
            value = sqrt(1 - M2.^2) ./ (beta .* (R.^2 .* M2.^2 + (1 - M2.^2) .* zeta.^2)) .* (-zeta + (beta .* R.^2 .* M2.^2) ./ sqrt(1 - M2.^2) .* sqrt( (1 - M2.^2) ./  (R.^2 .* M2.^2) .* zeta.^2 + ((R .* beta .* M2).^2 - 1) ./ (R .* beta .* M2).^2));
        end
        
        function value = zeta2theta(R, M2, beta, zeta)
            % Definitions
            costheta = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.costheta(R, M2, beta, zeta);

            % Compute theta
            value = acos(costheta);
        end
        
        function [thetaCritical1, thetaCritical2] = thetaCritical(R, M2, beta)
            thetaCritical1 = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.zeta2theta(R, M2, beta, -1);
            thetaCritical2 = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.zeta2theta(R, M2, beta, 1);
        end
        
        function value = theta2zeta(R, M2, beta, theta)
            % Definitions
            M1 = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.M1(R, M2, beta);

            % Compute zeta
            value = (cos(theta) - 1 ./ M1) .* R .* M2 ./ sqrt(1 - M2.^2) ./ sin(theta);
        end
        
        function value = sigma_a(R, M2, Gammas)
            value =  (R  ./ (R  - 1)) .* ((1 - Gammas) ./ (2*M2));
        end
        
        function value = sigma_b(M2, Gammas)
            value =  (1 + Gammas) ./ (2*M2);
        end
        
        function value = sigma_c(R, M2, Gammas)
            % Definitions
            sigma_a = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.sigma_a(R, M2, Gammas);

            % Compute sigma_c
            value =  (((M2.^2 .* (R - 1))) ./ (1 - M2.^2)) .* sigma_a;
        end
        
        function value = sigma_d(R, M2, Gammas, Gammas1)
            value =  M2 .* R .* (Gammas + Gammas1) ./ (2 * Gammas1);
        end
        
        function value = sigma_e(R, Gammas, Gammas3, beta)
            value =  R .* Gammas.^2 ./ Gammas3 - 1 ./ beta.^2;
        end
        
        function value = alpha1(R, M2, Gammas, Gammas1, Gammas3, beta, zeta)
            % Definitions
            sigma_d = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.sigma_d(R, M2, Gammas, Gammas1);
            sigma_e = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.sigma_e(R, Gammas, Gammas3, beta);
            costheta = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.costheta(R, M2, beta, zeta);
            term1 = ( 2 + beta.^2 .* (2 .* M2 .* R .* sigma_d + sigma_e) ) ./ (M2 .* (M2.^2 - 1) .* R .* beta.^2);
            term2 = 2 .* costheta ./ ((M2.^2 - 1) .* beta);
            
            % Compute alpha1
            value = (1 - M2.^2) .* (zeta.^2 / 2) .* (term1 - term2);
        end
        
        function value = alpha2(R, M2, Gammas, Gammas1, Gammas3, beta, zeta)
            % Definitions
            sigma_d = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.sigma_d(R, M2, Gammas, Gammas1);
            sigma_e = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.sigma_e(R, Gammas, Gammas3, beta);
            costheta = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.costheta(R, M2, beta, zeta);
            sintheta = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.sintheta(R, M2, beta, zeta);
            term1 = M2 .* (-2 .* M2.^2 .* R + (2 .* M2 .* R .* (Gammas1 - Gammas) .* sigma_d) ./ (Gammas1 + Gammas) - sigma_e);
            term2 = (2 .* M2.^2 .* (R - 1) .* costheta) ./ beta;
            term3 = (M2 .* zeta .* sintheta) ./ (sqrt(1 - M2.^2) .* beta);
            
            % Compute alpha2
            value = (1 ./ (2 .* (1 - M2.^2))) .* (term1 - term2) + term3;
        end
        
        function value = alpha(R, M2, Gammas, Gammas1, Gammas3, beta, zeta)
            % Definitions
            alpha1 = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.alpha1(R, M2, Gammas, Gammas1, Gammas3, beta, zeta);
            alpha2 = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.alpha2(R, M2, Gammas, Gammas1, Gammas3, beta, zeta);

            % Compute alpha
            value = alpha1 + alpha2;
        end
        
        function value = pi_l1(R, M2, Gammas, Gammas1, Gammas3, beta, zeta)
            % Definitions
            sigma_b = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.sigma_b(M2, Gammas);
            sigma_c = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.sigma_c(R, M2, Gammas);
            alpha = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.alpha(R, M2, Gammas, Gammas1, Gammas3, beta, zeta);
            
            % Compute the denominator
            denom = zeta.^2 .* (1 - zeta.^2) + (sigma_b .* zeta.^2 - sigma_c).^2;
            
            % Compute the numerator
            num = -alpha .* (sigma_b .* zeta.^2 - sigma_c);
            
            % Compute the value
            value = num ./ denom;
        end
        
        function value = pi_l2(R, M2, Gammas, Gammas1, Gammas3, beta, zeta)
            % Definitions
            sigma_b = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.sigma_b(M2, Gammas);
            sigma_c = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.sigma_c(R, M2, Gammas);
            alpha = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.alpha(R, M2, Gammas, Gammas1, Gammas3, beta, zeta);

            % Compute the denominator
            denom = zeta.^2 .* (1 - zeta.^2) + (sigma_b .* zeta.^2 - sigma_c).^2;
            
            % Compute the numerator
            num = (alpha) .* (zeta .* sqrt(1 - zeta.^2));
            
            % Compute the value
            value = num ./ denom;
        end
        
        function value = pi_s(R, M2, Gammas, Gammas1, Gammas3, beta, zeta)
            % Definitions
            sigma_b = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.sigma_b(M2, Gammas);
            sigma_c = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.sigma_c(R, M2, Gammas);
            alpha = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.alpha(R, M2, Gammas, Gammas1, Gammas3, beta, zeta);

            % Compute the denominator
            denom = zeta .* sqrt(zeta.^2 - 1) + (sigma_b .* zeta.^2 - sigma_c);
            
            % Compute the numerator
            num = -alpha;
            
            % Compute the value
            value = num ./ denom;
        end
        
        function value = Delta_pi(R, M2, Gammas, Gammas1, Gammas3, beta, zeta)

            if zeta > 1
                % Compute the value for shortwave (zeta > 1)
                value = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.pi_s(R, M2, Gammas, Gammas1, Gammas3, beta, zeta);
                return
            end
        
            % Definitions
            pi_l1 = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.pi_l1(R, M2, Gammas, Gammas1, Gammas3, beta, zeta);
            pi_l2 = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.pi_l2(R, M2, Gammas, Gammas1, Gammas3, beta, zeta);

            % Compute the value for longwave (zeta < 1)
            value = sqrt(pi_l1.^2 + pi_l2.^2);
        end
        
        function value = ka(M2, zeta)
            value = (zeta .* M2 - sqrt(zeta.^2 - 1)) ./ sqrt(1 - M2.^2);
        end
        
        function value = wa(M2, zeta)
            value = (zeta - M2 .* sqrt(zeta.^2 - 1)) ./ sqrt(1 - M2.^2);
        end
        
        function value = ws(M2, zeta)
            % Definitions
            wa = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.wa(M2, zeta);
            ka = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.ka(M2, zeta);

            % Compute ws
            value = wa - ka .* M2;
        end
        
        function value = Xs(R, M2, Gammas, Gammas1, Gammas3, beta, zeta)
            % Definitions
            sigma_c = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.sigma_c(R, M2, Gammas);
            sigma_d = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.sigma_d(R, M2, Gammas, Gammas1);
            pi_s = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.pi_s(R, M2, Gammas, Gammas1, Gammas3, beta, zeta);
            ws = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.ws(M2, zeta);

            % Compute Xs
            value = ((1 - M2.^2) .* sigma_c) ./ (M2.^2 .* (R - 1)) .* pi_s ./ ws ...
                    - R ./ (ws .* (R - 1)) .* ((Gammas1 - Gammas) ./ (Gammas1 + Gammas) .* sigma_d - M2) ...
                    + 1 ./ (beta .* ws);
        end
        
        function value = Xl1(R, M2, Gammas, Gammas1, Gammas3, beta, zeta)
            % Definitions
            sigma_c = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.sigma_c(R, M2, Gammas);
            sigma_d = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.sigma_d(R, M2, Gammas, Gammas1);
            pi_l1 = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.pi_l1(R, M2, Gammas, Gammas1, Gammas3, beta, zeta);
            ws = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.ws(M2, zeta);

            % Compute Xl1
            value = ((1 - M2.^2) .* sigma_c) ./ (M2.^2 .* (R - 1)) .* pi_l1 ./ ws ...
                - R ./ (ws .* (R - 1)) .* ((Gammas1 - Gammas) ./ (Gammas1 + Gammas) .* sigma_d - M2) ...
                + 1 ./ (beta .* ws);
        end
        
        function value = Xl2(R, M2, Gammas, Gammas1, Gammas3, beta, zeta)
            % Definitions
            sigma_c = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.sigma_c(R, M2, Gammas);
            pi_l2 = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.pi_l2(R, M2, Gammas, Gammas1, Gammas3, beta, zeta);
            ws = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.ws(M2, zeta);

            % Compute Xl2
            value = ((1 - M2.^2) .* sigma_c) ./ (M2.^2 .* (R - 1)) .* pi_l2 ./ ws;
        end
        
        function value = Xl(R, M2, Gammas, Gammas1, Gammas3, beta, zeta)
            % Definitions
            Xl1 = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.Xl1(R, M2, Gammas, Gammas1, Gammas3, beta, zeta);
            Xl2 = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.Xl2(R, M2, Gammas, Gammas1, Gammas3, beta, zeta);

            % Compute Xl
            value = sqrt(Xl1.^2 + Xl2.^2);
        end
        
        function value = Delta_ua(R, M2, Gammas, Gammas1, Gammas3, beta, zeta)
            % Definitions
            ka = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.ka(M2, zeta);
            wa = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.wa(M2, zeta);
            pi_s = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.pi_s(R, M2, Gammas, Gammas1, Gammas3, beta, zeta);

            % Compute Delta_ua
            value = (ka ./ wa) .* pi_s;
        end
        
        function value = Delta_va(R, M2, Gammas, Gammas1, Gammas3, beta, zeta)
            % Definitions
            wa = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.wa(M2, zeta);
            pi_s = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.pi_s(R, M2, Gammas, Gammas1, Gammas3, beta, zeta);

            % Compute Delta_va

            value = (1 ./ wa) .* pi_s;
        end
        
        function value = Delta_rho_l1(R, M2, Gammas, Gammas1, Gammas3, beta, zeta)
            % Definitions
            pi_l1 = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.pi_l1(R, M2, Gammas, Gammas1, Gammas3, beta, zeta);

            % Compute Delta_rho_l1
            value = (Gammas ./ M2.^2 - 1) .* pi_l1 - R .* Gammas ./ Gammas1;
        end
        
        function value = Delta_rho_l2(R, M2, Gammas, Gammas1, Gammas3, beta, zeta)
            % Definitions
            pi_l2 = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.pi_l2(R, M2, Gammas, Gammas1, Gammas3, beta, zeta);

            % Compute Delta_rho_l2
            value = (Gammas ./ M2.^2 - 1) .* pi_l2;
        end
        
        function value = Delta_rho_s(R, M2, Gammas, Gammas1, Gammas3, beta, zeta)
            % Definitions
            pi_s = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.pi_s(R, M2, Gammas, Gammas1, Gammas3, beta, zeta);

            % Compute Delta_rho_s
            value = (Gammas ./ M2.^2 - 1) .* pi_s - R .* Gammas ./ Gammas1;
        end
        
        function value = Delta_rho(R, M2, Gammas, Gammas1, Gammas3, beta, zeta)

            if zeta > 1
                % Compute the value for shortwave (zeta > 1)
                value = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.Delta_rho_s(R, M2, Gammas, Gammas1, Gammas3, beta, zeta);
                return
            end
        
            % Definitions
            Delta_rho_l1 = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.Delta_rho_l1(R, M2, Gammas, Gammas1, Gammas3, beta, zeta);
            Delta_rho_l2 = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.Delta_rho_l2(R, M2, Gammas, Gammas1, Gammas3, beta, zeta);

            % Compute Delta_rho
            value = sqrt(Delta_rho_l1.^2 + Delta_rho_l2.^2);
        end
        
        function value = Omega_1(R, M2, beta, zeta)
            value = R ./ beta .* (1 + zeta.^2 .* ((1 - M2.^2) ./ (R.^2 .* M2.^2)));
        end
        
        function value = Omega_2(R, M2, Gammas)
            value = ((R - 1) .* (1 - Gammas)) ./ (2 .* M2);
        end
        
        function value = Omega_3(R, M2, Gammas, Gammas1)
            % Definitions
            sigma_d = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.sigma_d(R, M2, Gammas, Gammas1);

            % Compute Omega_3
            value = -(R - 1) .* (R .* M2 - sigma_d);
        end
        
        function value = Omega_4(R, M2, Gammas)
            value = - combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.Omega_2(R, M2, Gammas);
        end
        
        function value = Delta_Omega_l1(R, M2, Gammas, Gammas1, Gammas3, beta, zeta)
            % Definitions
            Omega_1 = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.Omega_1(R, M2, beta, zeta);
            Omega_2 = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.Omega_2(R, M2, Gammas);
            Omega_3 = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.Omega_3(R, M2, Gammas, Gammas1);
            pi_l1 = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.pi_l1(R, M2, Gammas, Gammas1, Gammas3, beta, zeta);

            % Compute Delta_Omega_l1
            value = Omega_1 + Omega_2 .* pi_l1 + Omega_3;
        end
        
        function value = Delta_Omega_l2(R, M2, Gammas, Gammas1, Gammas3, beta, zeta)
            % Definitions
            Omega_2 = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.Omega_2(R, M2, Gammas);
            pi_l2 = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.pi_l2(R, M2, Gammas, Gammas1, Gammas3, beta, zeta);

            % Compute Delta_Omega_l2
            value = Omega_2 .* pi_l2;
        end
        
        function value = Delta_Omega_s(R, M2, Gammas, Gammas1, Gammas3, beta, zeta)
            % Definitions
            Omega_1 = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.Omega_1(R, M2, beta, zeta);
            Omega_2 = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.Omega_2(R, M2, Gammas);
            Omega_3 = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.Omega_3(R, M2, Gammas, Gammas1);
            pi_s = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.pi_s(R, M2, Gammas, Gammas1, Gammas3, beta, zeta);

            % Compute Delta_Omega_s
            value = Omega_1 + Omega_2 .* pi_s + Omega_3;
        end
        
        function value = Delta_Omega_l(R, M2, Gammas, Gammas1, Gammas3, beta, zeta)
            % Definitions
            Delta_Omega_l1 = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.Delta_Omega_l1(R, M2, Gammas, Gammas1, Gammas3, beta, zeta);
            Delta_Omega_l2 = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.Delta_Omega_l2(R, M2, Gammas, Gammas1, Gammas3, beta, zeta);

            % Compute Delta_Omega_l
            value = sqrt(Delta_Omega_l1.^2 + Delta_Omega_l2.^2);
        end
        
        function value = Delta(M2, zeta)
            value = 1 + ((1 - M2.^2) ./ M2.^2) .* zeta.^2;
        end
        
        function value = Delta_u_l1(R, M2, Gammas, Gammas1, Gammas3, beta, zeta)
            % Definitions
            Delta_Omega_l1 = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.Delta_Omega_l1(R, M2, Gammas, Gammas1, Gammas3, beta, zeta);
            Delta = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.Delta(M2, zeta);

            % Compute Delta_u_l1
            value = Delta_Omega_l1 ./ Delta;
        end
        
        function value = Delta_u_l2(R, M2, Gammas, Gammas1, Gammas3, beta, zeta)
            % Definitions
            Delta_Omega_l2 = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.Delta_Omega_l2(R, M2, Gammas, Gammas1, Gammas3, beta, zeta);
            Delta = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.Delta(M2, zeta);

            % Compute Delta_u_l2
            value = Delta_Omega_l2 ./ Delta;
        end
        
        function value = Delta_u_s(R, M2, Gammas, Gammas1, Gammas3, beta, zeta)
            % Definitions
            Delta_Omega_s = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.Delta_Omega_s(R, M2, Gammas, Gammas1, Gammas3, beta, zeta);
            Delta = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.Delta(M2, zeta);

            % Compute Delta_u_s
            value = Delta_Omega_s ./ Delta;
        end
        
        function value = Delta_u(R, M2, Gammas, Gammas1, Gammas3, beta, zeta)

            if zeta > 1
                % Compute the value for shortwave (zeta > 1)
                value = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.Delta_u_s(R, M2, Gammas, Gammas1, Gammas3, beta, zeta);
                return
            end
            
            % Definitions
            Delta_u_l1 = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.Delta_u_l1(R, M2, Gammas, Gammas1, Gammas3, beta, zeta);
            Delta_u_l2 = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.Delta_u_l2(R, M2, Gammas, Gammas1, Gammas3, beta, zeta);

            % Compute Delta_u
            value = sqrt(Delta_u_l1.^2 + Delta_u_l2.^2);
        end
        
        function value = Delta_v_l1(R, M2, Gammas, Gammas1, Gammas3, beta, zeta)
            % Definitions
            Delta_Omega_l1 = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.Delta_Omega_l1(R, M2, Gammas, Gammas1, Gammas3, beta, zeta);
            Delta = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.Delta(M2, zeta);

            % Compute Delta_v_l1
            value = zeta.* ((sqrt(1 - M2.^2)) ./ M2) .* (Delta_Omega_l1 ./ Delta);
        end
        
        function value = Delta_v_l2(R, M2, Gammas, Gammas1, Gammas3, beta, zeta)
            % Definitions
            Delta_Omega_l2 = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.Delta_Omega_l2(R, M2, Gammas, Gammas1, Gammas3, beta, zeta);
            Delta = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.Delta(M2, zeta);

            % Compute Delta_v_l2
            value = zeta.* ((sqrt(1 - M2.^2)) ./ M2) .* (Delta_Omega_l2 ./ Delta);
        end
        
        function value = Delta_v_s(R, M2, Gammas, Gammas1, Gammas3, beta, zeta)
            % Definitions
            Delta_Omega_s = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.Delta_Omega_s(R, M2, Gammas, Gammas1, Gammas3, beta, zeta);
            Delta = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.Delta(M2, zeta);

            % Compute Delta_v_s
            value = zeta.* ((sqrt(1 - M2.^2)) ./ M2) .* (Delta_Omega_s ./ Delta);
        end
        
        function value = Delta_v(R, M2, Gammas, Gammas1, Gammas3, beta, zeta)

            if zeta > 1
                % Compute the value for shortwave (zeta > 1)
                value = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.Delta_v_s(R, M2, Gammas, Gammas1, Gammas3, beta, zeta);
                return
            end
            
            % Definitions
            Delta_v_l1 = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.Delta_v_l1(R, M2, Gammas, Gammas1, Gammas3, beta, zeta);
            Delta_v_l2 = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.Delta_v_l2(R, M2, Gammas, Gammas1, Gammas3, beta, zeta);

            % Compute Delta_v
            value = sqrt(Delta_v_l1.^2 + Delta_v_l2.^2);
        end
        
        function value = b01(R, M2, Gammas, Gammas1, Gammas3, beta, zeta)
            % Definitions
            sigma_b = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.sigma_b(M2, Gammas);
            sigma_c = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.sigma_c(R, M2, Gammas);
            alpha = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.alpha(R, M2, Gammas, Gammas1, Gammas3, beta, zeta);

            % Compute the numerator (sum of alpha1 and alpha2)
            num = alpha .* (-sigma_b .* zeta.^2 + sigma_c);
            
            % Compute the denominator
            denom = zeta.^2 .* (1 - zeta.^2) + (-sigma_b .* zeta.^2 + sigma_c).^2;
            
            % Compute the value
            value = num ./ denom;
        end

        function value = b02(R, M2, Gammas, Gammas1, Gammas3, beta, zeta)
            % Definitions
            sigma_b = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.sigma_b(M2, Gammas);
            sigma_c = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.sigma_c(R, M2, Gammas);
            alpha = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.alpha(R, M2, Gammas, Gammas1, Gammas3, beta, zeta);

            % Compute the numerator (sum of alpha1 and alpha2)
            num = alpha .* zeta .* sqrt(1 - zeta.^2);
            
            % Compute the denominator
            denom = zeta.^2 .* (1 - zeta.^2) + (-sigma_b .* zeta.^2 + sigma_c).^2;
            
            % Compute the value
            value = num ./ denom;
        end

        function value = pe0(R, M2, Gammas, Gammas1, Gammas3, beta, zeta)
            % Definitions
            sigma_b = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.sigma_b(M2, Gammas);
            sigma_c = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.sigma_c(R, M2, Gammas);
            alpha = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.alpha(R, M2, Gammas, Gammas1, Gammas3, beta, zeta);

            % Compute the numerator (sum of alpha1 and alpha2)
            num = -alpha;
            
            % Compute the denominator
            denom = zeta .* sqrt(zeta.^2 - 1) + sigma_b .* zeta.^2 - sigma_c;
            
            % Compute the value
            value = num ./ denom;
        end

        function value = e0(R, M2, Gammas, Gammas1, Gammas3, beta, zeta)
            % Definitions
            sigma_b = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.sigma_b(M2, Gammas);
            sigma_c = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.sigma_c(R, M2, Gammas);
            alpha = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.alpha(R, M2, Gammas, Gammas1, Gammas3, beta, zeta);

            % Compute the numerator (sum of alpha1 and alpha2)
            num = -alpha;
            
            % Compute the denominator
            denom = abs(zeta) .* sqrt(zeta.^2 - 1) + sigma_b .* zeta.^2 - sigma_c;
            
            % Compute the value
            value = num ./ denom;
        end
        
        function value = ne0(R, M2, Gammas, Gammas1, Gammas3, beta, zeta)
            % Definitions
            sigma_b = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.sigma_b(M2, Gammas);
            sigma_c = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.sigma_c(R, M2, Gammas);
            alpha = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.alpha(R, M2, Gammas, Gammas1, Gammas3, beta, zeta);

            % Compute the numerator (sum of alpha1 and alpha2)
            num = -alpha;
            
            % Compute the denominator
            denom = -zeta .* sqrt(zeta.^2 - 1) + sigma_b .* zeta.^2 - sigma_c;
            
            % Compute the value
            value = num ./ denom;
        end
        
        function value = Qsrot(R, M2, Gammas, Gammas1, Gammas3, beta, zeta)
            % Definitions
            Omega_2 = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.Omega_2(R, M2, Gammas);
            e0 = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.e0(R, M2, Gammas, Gammas1, Gammas3, beta, zeta);
            Delta = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.Delta(M2, zeta);

            % Compute Qsrot
            value = (Omega_2 .* (e0 - 1)) ./ Delta;
        end
        
        function value = Qlrot(R, M2, Gammas, Gammas1, Gammas3, beta, zeta)
            % Definitions
            Omega_2 = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.Omega_2(R, M2, Gammas);
            Delta = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.Delta(M2, zeta);
            pi_l1 = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.pi_l1(R, M2, Gammas, Gammas1, Gammas3, beta, zeta);
            pi_l2 = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.pi_l2(R, M2, Gammas, Gammas1, Gammas3, beta, zeta);

            % Compute Qlrot
            value = -(Omega_2 ./ Delta) .* sqrt((pi_l1 - 1).^2 + pi_l2.^2);
        end
        
        function value = Qax(R, M2, Gammas, Gammas1, Gammas3, beta, zeta)
            % Definitions
            e0 = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic.e0(R, M2, Gammas, Gammas1, Gammas3, beta, zeta);

            % Compute Qax
            value = ( (sqrt(zeta.^2 - 1) - abs(zeta) .* M2) ./ (M2 .* sqrt(zeta.^2 - 1) - abs(zeta)) ) .* e0;
        end

    end

end
