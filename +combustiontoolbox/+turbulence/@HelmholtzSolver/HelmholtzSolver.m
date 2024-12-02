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
    % Args:
    %     file_location (char): Path to the data .hdf file
    %     file_location_nodes (char): Path to the grid .hdf file
    %     T_ref (float): Temperature of reference [K]
    %     mu_ref (float): Dynamic viscosity of reference [kg/(m-s)] or [Pa-s]
    %
    % Example:
    %     solver = HelmholtzSolver();
    %
    % References:
    %   [1] Johnson, S. G. (2011). Notes on FFT-based differentiation.
    %       MIT Applied Mathematics, Tech. Rep. 
    %       Available: http://math.mit.edu/~stevenj/fft-deriv.pdf
    %   [2] Xun Shi, Helmholtz-Hodge decomposition using fft (Python),
    %       Available: https://github.com/shixun22/helmholtz
    %
    %
    % @author: Alberto Cuadra Lara
    %          Postdoctoral researcher - Group Fluid Mechanics
    %          Universidad Carlos III de Madrid

    properties
        tol0 = 1e-3             % Tolerance for checks
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
            addParameter(p, 'FLAG_CHECKS', obj.FLAG_CHECKS, @(x) islogical(x));
            addParameter(p, 'FLAG_TIME', obj.FLAG_TIME, @(x) islogical(x));
            addParameter(p, 'FLAG_REPORT', obj.FLAG_REPORT, @(x) islogical(x));
            addParameter(p, 'plotConfig', defaultPlotConfig, @(x) isa(x, 'combustiontoolbox.utils.display.PlotConfig'));
            parse(p, varargin{:});

            % Set properties
            obj.tol0 = p.Results.tol0;
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
            %     * rho (float): Density field
            %
            % Returns:
            %     solenoidal (VelocityField): Struct with fields (u, v, w) containing the solenoidal velocity components (fluctuations)
            %     dilatational (VelocityField): Struct with fields (u, v, w) containing the dilatational velocity components (fluctuations)
            %     velocity (VelocityField): Struct with fields (u, v, w) containing the velocity components (fluctuations)
            %     STOP (float): Relative error doing the decomposition
            %
            % Examples:
            %     * [solenoidal, dilatational, velocity, STOP] = solve(obj, velocity)
            %     * [solenoidal, dilatational, velocity, STOP] = solve(obj, velocity, 'rho', rho)
            
            % Import packages
            import combustiontoolbox.common.Units.convertData2VelocityField

            % Timer
            obj.time = tic;

            % Parse input arguments
            p = inputParser;
            addOptional(p, 'rho', [], @(x) isnumeric(x));
            parse(p, varargin{:});

            % Set properties
            rho = p.Results.rho;
            
            % Reshape velocity input and compute fluctuations
            velocity = convertData2VelocityField(velocity);
            velocity = obj.computeFluctuations(velocity, rho);

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
            U = fftn(velocity.u);
            V = fftn(velocity.v);
            W = fftn(velocity.w);

            % Compute wave numbers
            [KX, KY, KZ] = obj.computeWaveNumbers(size(velocity.u));

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
            [KX, KY, KZ] = obj.computeWaveNumbers(size(velocity.u));
            
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
                fprintf('Error!\nReconstructed field does not match original field: %.2e\n', diffMax);
                STOP = diffMax;
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
    
        function velocity = computeFluctuations(velocity, rho)
            % Compute fluctuating velocity components
            %
            % Args:
            %     velocity (VelocityField): VelocityField instance with fields (u, v, w) containing the velocity components
            %     rho (float): Density field
            %
            % Returns:
            %     velocity (struct): Struct with fields (u, v, w) containing the fluctuating velocity components (fluctuations)
            
            % Shallow copy of the velocity field
            velocity = velocity.copy();

            % For compressible flows
            if ~isempty(rho)
                rhoMean = mean(rho, 'all');
                rhou = mean(rho .* velocity.u, 'all') / rhoMean;
                rhov = mean(rho .* velocity.v, 'all') / rhoMean;
                rhow = mean(rho .* velocity.w, 'all') / rhoMean;

                velocity.u = sqrt(rho) .* (velocity.u - rhou);
                velocity.v = sqrt(rho) .* (velocity.v - rhov);
                velocity.w = sqrt(rho) .* (velocity.w - rhow);
                return
            end

            % For incompressible flows
            velocity.u = velocity.u - mean(velocity.u, 'all');
            velocity.v = velocity.v - mean(velocity.v, 'all');
            velocity.w = velocity.w - mean(velocity.w, 'all');
        end

        function [KX, KY, KZ] = computeWaveNumbers(sz)
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
            %     [KX, KY, KZ] = computeWaveNumbers(sz)

            % Import packages
            import combustiontoolbox.utils.math.fftfreq

            kx = fftfreq(sz(1));
            ky = fftfreq(sz(2));
            kz = fftfreq(sz(3));
            [KX, KY, KZ] = ndgrid(kx, ky, kz);
        end
        
    end

end