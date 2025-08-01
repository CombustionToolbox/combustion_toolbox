classdef HelmholtzSolver < handle
    % The :mat:func:`HelmholtzSolver` class computes the Helmholtz-Hodge 
    % decomposition of a 3D velocity field into its solenoidal and
    % dilatational parts using fast Fourier transform [1].
    % 
    % The decomposition is performed with the spectral method, which is
    % only suitable for relatively smooth fields, i.e., with little power
    % on small scales. The code assumes that the grid is uniform with
    % dx = dy = dz.
    %
    % This code is based on Ref. [2] and has been rewritten in MATLAB 
    % with some modifications.
    %
    % Notes:
    %     For even NX, NY, and NZ, decomposed fields can be complex,
    %     with the imaginary part coming from the real part of the kmode
    %     at Nyquist frequency. In principle, the Nyquist frequency
    %     kmode should be dropped when doing the first derivatives to
    %     maintain symmetry. See footnote on page 4 of [2]. However,
    %     when the field is smooth enough, the imaginary part caused by
    %     the Nyquist frequency kmode should be negligible.
    %
    % References:
    %   [1] Johnson, S. G. (2011). Notes on FFT-based differentiation.
    %       MIT Applied Mathematics, Tech. Rep. 
    %       Available: http://math.mit.edu/~stevenj/fft-deriv.pdf
    %   [2] Xun Shi, Helmholtz-Hodge decomposition using fft (Python),
    %       Available: https://github.com/shixun22/helmholtz

    properties
        tol0 = 1e-3             % Tolerance for checks
        FLAG_WEIGHTED = true    % Flag to compute weighted velocity field as u * sqrt(rho)
        FLAG_CHECKS = true      % Flag to perform checks
        FLAG_TIME = true        % Flag to print elapsed time
        FLAG_REPORT = false     % Flag to print predefined plots
        time                    % Elapsed time
        plotConfig              % PlotConfig object
    end

    methods 

        function obj = HelmholtzSolver(varargin)
            % Constructor 
            defaultPlotConfig = combustiontoolbox.utils.display.PlotConfig();

            % Parse input arguments
            p = inputParser;
            addParameter(p, 'tol0', obj.tol0, @(x) isnumeric(x));
            addParameter(p, 'FLAG_WEIGHTED', obj.FLAG_WEIGHTED, @(x) islogical(x));
            addParameter(p, 'FLAG_CHECKS', obj.FLAG_CHECKS, @(x) islogical(x));
            addParameter(p, 'FLAG_TIME', obj.FLAG_TIME, @(x) islogical(x));
            addParameter(p, 'FLAG_REPORT', obj.FLAG_REPORT, @(x) islogical(x));
            addParameter(p, 'plotConfig', defaultPlotConfig, @(x) isa(x, 'combustiontoolbox.utils.display.PlotConfig'));
            parse(p, varargin{:});

            % Set properties
            obj.tol0 = p.Results.tol0;
            obj.FLAG_WEIGHTED = p.Results.FLAG_WEIGHTED;
            obj.FLAG_CHECKS = p.Results.FLAG_CHECKS;
            obj.FLAG_TIME = p.Results.FLAG_TIME;
            obj.FLAG_REPORT = p.Results.FLAG_REPORT;
            obj.plotConfig = p.Results.plotConfig;
        end

        function obj = set(obj, property, value, varargin)
            % Set properties of the HelmholtzSolver object
            %
            % Args:
            %     obj (HelmholtzSolver): HelmholtzSolver object
            %     property (char): Property name
            %     value (float): Property value
            %
            % Optional Args:
            %     * property (char): Property name
            %     * value (float): Property value
            %
            % Returns:
            %     obj (HelmholtzSolver): HelmholtzSolver object with updated properties
            %
            % Examples:
            %     * set(HelmholtzSolver(), 'tol0', 1e-10)
            %     * set(HelmholtzSolver(), 'tol0', 1e-10, 'FLAG_CHECKS', false)
            
            varargin = [{property, value}, varargin{:}];

            for i = 1:2:length(varargin)
                % Assert that the property exists
                assert(isprop(obj, varargin{i}), 'Property not found');

                % Set property
                obj.(varargin{i}) = varargin{i + 1};
            end

        end

        function [solenoidal, dilatational, velocity, STOP] = solve(obj, velocity, varargin)
            % Compute Helmholtz decomposition of the velocity field
            %
            % Args:
            %     obj (HelmholtzSolver): HelmholtzSolver object
            %     velocity (VelocityField): Velocity field as a VelocityField object, struct, or 4D matrix
            %
            % Optional Args:
            %     * density (float): Density field
            %
            % Returns:
            %     solenoidal (VelocityField): VelocityField object with fields (u, v, w) containing the solenoidal velocity components (fluctuations)
            %     dilatational (VelocityField): VelocityField with fields (u, v, w) containing the dilatational velocity components (fluctuations)
            %     velocity (VelocityField): VelocityField with fields (u, v, w) containing the velocity components (fluctuations)
            %     STOP (float): Relative error doing the decomposition
            %
            % Shortcuts:
            %     * [solenoidal, dilatational, velocity, STOP] = solve(obj, velocity, rho)
            %
            % Examples:
            %     * [solenoidal, dilatational, velocity, STOP] = solve(obj, velocity)
            %     * [solenoidal, dilatational, velocity, STOP] = solve(obj, velocity, rho)
            %     * [solenoidal, dilatational, velocity, STOP] = solve(obj, velocity, 'density', rho)
            
            % Import packages
            import combustiontoolbox.common.Units.convertData2VelocityField

            % Timer
            obj.time = tic;
            
            % Shortcut: solve(obj, velocity, 'x')
            if isscalar(varargin)
                varargin = {'density', varargin{1}};
            end

            % Parse input arguments
            p = inputParser;
            addOptional(p, 'density', [], @(x) isnumeric(x));
            parse(p, varargin{:});

            % Set properties
            rho = p.Results.density;
            
            % Reshape velocity input and compute fluctuations
            velocity = convertData2VelocityField(velocity);
            velocity = getFluctuations(velocity, 'density', rho, 'weighted', obj.FLAG_WEIGHTED);

            % Solve the Helmholtz equation
            [solenoidal, dilatational] = obj.decomposition(velocity);
            
            % Check results
            STOP = obj.check(velocity, solenoidal, dilatational);

            % Time elapsed
            obj.time = toc(obj.time);
            
            % Print elapsed time
            printTime(obj);
        end

        function [solenoidal, dilatational] = decomposition(obj, velocity)
            % Compute Helmholtz decomposition of the velocity field (fluctuations)
            %
            % Args:
            %     velocity (VelocityField): VelocityField instance with fields (u, v, w) containing the velocity components (fluctuations)
            %
            % Returns:
            %     solenoidal (VelocityField): VelocityField instance with fields (u, v, w) containing the solenoidal velocity components (fluctuations)
            %     dilatational (VelocityField): VelocityField instance with fields (u, v, w) containing the dilatational velocity components (fluctuations)
            
            % Import packages
            import combustiontoolbox.turbulence.VelocityField

            % Get N-D fast Fourier transform (fft)
            [U, V, W] = getFFT(velocity);

            % Compute wave numbers
            [KX, KY, KZ] = obj.getWaveNumbers(size(velocity.u));

            % Compute k^2, avoiding division by zero
            K2 = KX.^2 + KY.^2 + KZ.^2;
            K2(K2 == 0) = 1;

            % Compute velocity divergence
            div = (U .* KX + V .* KY + W .* KZ);
            clear U V W % Free memory

            % Compute the Helmholtz decomposition (dilatational)
            H = div ./ K2;
            clear div K2 % Free memory

            % Compute dilatational components
            dilatational = VelocityField(real(ifftn(H .* KX)), ...
                                         real(ifftn(H .* KY)), ...
                                         real(ifftn(H .* KZ)));
            clear H KX KY KZ; % Free memory

            % Compute solenoidal components
            solenoidal = VelocityField(velocity.u - dilatational.u, ...
                                       velocity.v - dilatational.v, ...
                                       velocity.w - dilatational.w);
        end

        function [STOP, status] = check(obj, velocity, solenoidal, dilatational)
            % Check the results of the Helmholtz decomposition
            %
            % Args:
            %     obj (HelmholtzSolver): HelmholtzSolver object
            %     velocity (VelocityField): VelocityField instance with fields (u, v, w) containing the velocity components (fluctuations)
            %     solenoidal (VelocityField): VelocityField instance with fields (u, v, w) containing the solenoidal velocity components (fluctuations)
            %     dilatational (VelocityField): VelocityField instance with fields (u, v, w) containing the dilatational velocity components (fluctuations)
            %
            % Returns:
            %     STOP (float): Relative error doing the decomposition
            %     status (bool): Flag indicating if the checks passed

            if ~obj.FLAG_CHECKS
                status = [];
                STOP = [];
                return
            end

            % Validate Helmholtz decomposition results
            fprintf('Performing checks... ');
            
            % Compute wave numbers
            [KX, KY, KZ] = obj.getWaveNumbers(size(velocity.u));
            
            % Compute magnitude of the original velocity field for relative error
            velocityMagnitude = globalMagnitude(velocity);

            % Check if the solenoidal part is divergence-free
            divSolenoidal = ifftn((fftn(solenoidal.u) .* KX + ...
                                   fftn(solenoidal.v) .* KY + ...
                                   fftn(solenoidal.w) .* KZ));
                                   
            divError = max(abs(divSolenoidal(:))) / velocityMagnitude;

            if divError > obj.tol0
                fprintf('Error!\nSolenoidal part divergence too high: %.2e\n', divError);
                STOP = divError;
                return;
            end

            % Check if the dilatational part is curl-free
            curlDilatational = ifftn((fftn(dilatational.w) .* KY - fftn(dilatational.v) .* KZ + ...
                                      fftn(dilatational.u) .* KZ - fftn(dilatational.w) .* KX + ...
                                      fftn(dilatational.v) .* KX - fftn(dilatational.u) .* KY) * 1i * 2 * pi);

            curlError = max(abs(curlDilatational(:))) / velocityMagnitude;
            
            if curlError > obj.tol0
                fprintf('Error!\nDilatational part curl too high: %.2e\n', curlError);
                STOP = curlError;
                return;
            end

            % Check if the solenoidal and dilatational parts sum up to the original field
            uRecon = solenoidal.u + dilatational.u;
            vRecon = solenoidal.v + dilatational.v;
            wRecon = solenoidal.w + dilatational.w;
            
            diffError = sqrt(sum((uRecon(:) - velocity.u(:)).^2 + ...
                                 (vRecon(:) - velocity.v(:)).^2 + ...
                                 (wRecon(:) - velocity.w(:)).^2)) / velocityMagnitude;
            
            if diffError > obj.tol0
                fprintf('Error!\nReconstructed field does not match original field: %.2e\n', diffError);
                STOP = diffError;
                return;
            end

            % Set flag to pass and stop
            STOP = max([divError, curlError, diffError]);
            status = true;
            fprintf('OK!\n');
        end

        function printTime(obj)
            % Print execution time
            %
            % Args:
            %     obj (EquilibriumSolver): Object of the class EquilibriumSolver
            
            if ~obj.FLAG_TIME
                return
            end

            % Definitions
            operationName = 'Helmholtz decomposition';

            % Print elapsed time
            fprintf('\nElapsed time for %s: %.5f seconds\n', operationName, obj.time);
        end

    end

    methods (Static) 

        function [KX, KY, KZ] = getWaveNumbers(sz)
            % Compute wave number grids for FFT
            %
            % Args:
            %     sz (float): Size of the 3D array
            %
            % Returns:
            %     KX (float): 3D array with the wave number in the x-direction
            %     KY (float): 3D array with the wave number in the y-direction
            %     KZ (float): 3D array with the wave number in the z-direction
            %
            % Example:
            %     [KX, KY, KZ] = getWaveNumbers(sz)

            % Import packages
            import combustiontoolbox.utils.math.fftfreq

            kx = fftfreq(sz(1));
            ky = fftfreq(sz(2));
            kz = fftfreq(sz(3));
            [KX, KY, KZ] = ndgrid(kx, ky, kz);
        end

        function [chi, chiVariance] = getChi(solenoidal, rho, pressure, sound_mean)
            % Compute the correlation between the entropic density and the soleonidal part of the velocity field
            %
            % Args:
            %     solenoidal (VelocityField): VelocityField instance with fields (u, v, w) containing the solenoidal velocity components (fluctuations)
            %     rho (float): Density field
            %     pressure (float): Pressure field
            %     sound_mean (float): Mean speed of sound
            %
            % Returns:
            %     chi (float): Correlation between the entropic density and the soleonidal part of the velocity field
            %     chiVariance (float): Variance of the correlation
            %
            % Note: The components of solenoidal refer to the fluctuations of the density weighted velocity field, i.e.,
            %       u_weighted = u * sqrt(rho), v_weighted = v * sqrt(rho), w_weighted = w * sqrt(rho).
            %
            % Example:
            %     [chi, chiVariance] = getChi(obj, solenoidal, rho, pressure, sound_mean)

            % Definitions
            rho_mean = mean(rho, 'all');
            pressure_mean = mean(pressure, 'all');
            delta_rho = rho - rho_mean;
            delta_pressure = pressure - pressure_mean;
            delta_u_solenoidal = solenoidal.u ./ sqrt(rho);

            % Compute acoustic and entropic fluctuations of the density
            delta_rho_acoustic = delta_pressure ./ sound_mean.^2;
            delta_rho_entropic = delta_rho - delta_rho_acoustic;
            
            % Compute root mean square values
            delta_u_solenoidal_rms = sqrt(mean((delta_u_solenoidal).^2, 'all'));
            delta_rho_entropic_rms = sqrt(mean(delta_rho_entropic.^2, 'all'));

            % Compute chi
            chi = - mean(solenoidal.u .* delta_rho_entropic, 'all') /  mean(delta_u_solenoidal.^2, 'all') * sound_mean / rho_mean;
            chiVariance = (delta_rho_entropic_rms ./ delta_u_solenoidal_rms * sound_mean / rho_mean)^2;
        end

    end

end