classdef VelocityField < handle & matlab.mixin.Copyable
    % This :mat:func:`VelocityField` class provides methods for computing
    % properties of a velocity field.
    
    properties
        u      % x-component of the velocity field
        v      % y-component of the velocity field
        w      % z-component of the velocity field
        x      % x-component field
        y      % y-component field
        z      % z-component field
    end

    properties(Dependent)
        size   % Size of the velocity field
    end

    properties(Access = private)
        defaultLenght = 2 * pi % Default length of the domain in each direction
        gridChecked = false    % Flag to check if the grid is defined
    end

    methods
        
        function obj = VelocityField(u, v, w, varargin)
            % Constructor for VelocityField
            %
            % Args:
            %     u (float): x-component of velocity
            %     v (float): y-component of velocity
            %     w (float): z-component of velocity
            %
            % Returns:
            %     obj (VelocityField): VelocityField object
            
            % Definitions
            defaultX = [];
            defaultY = [];
            defaultZ = [];

            % Parse inputs
            p = inputParser;
            addRequired(p, 'u', @isnumeric);
            addRequired(p, 'v', @isnumeric);
            addRequired(p, 'w', @isnumeric);
            addOptional(p, 'x', defaultX, @isnumeric);
            addOptional(p, 'y', defaultY, @isnumeric);
            addOptional(p, 'z', defaultZ, @isnumeric);
            parse(p, u, v, w, varargin{:});

            % Assign properties
            obj.u = p.Results.u;
            obj.v = p.Results.v;
            obj.w = p.Results.w;
            obj.x = p.Results.x;
            obj.y = p.Results.y;
            obj.z = p.Results.z;
        end

        function sz = get.size(obj)
            % Get size of the velocity field
            %
            % Args:
            %     obj (VelocityField): VelocityField object
            %
            % Returns:
            %     sz (vector): Size of the velocity field
            
            sz = size(obj.u); % Assume all components (u, v, w) have the same size
        end

        function set.x(obj, value)
            % Set x-component field
            %
            % Args:
            %     obj (VelocityField): VelocityField object
            %     value (float): x-component field
            
            obj.x = value;
            obj.gridChecked = false;
        end

        function set.y(obj, value)
            % Set y-component field
            %
            % Args:
            %     obj (VelocityField): VelocityField object
            %     value (float): y-component field
            
            obj.y = value;
            obj.gridChecked = false;
        end

        function set.z(obj, value)
            % Set z-component field
            %
            % Args:
            %     obj (VelocityField): VelocityField object
            %     value (float): z-component field
            
            obj.z = value;
            obj.gridChecked = false;
        end

        function obj = plus(obj1, obj2)
            % Overload the plus operator to add the velocity field from two VelocityField objects
            %
            % Args:
            %     obj1 (VelocityField): VelocityField object
            %     obj2 (VelocityField): VelocityField object
            %
            % Returns:
            %     obj (VelocityField): VelocityField object 

            obj = combustiontoolbox.turbulence.VelocityField(...
                    obj1.u + obj2.u, ...
                    obj1.v + obj2.v, ...
                    obj1.w + obj2.w, ...
                    obj1.x, obj1.y, obj1.z ...
                );
        end
        
        function velocity = getFluctuations(obj, varargin)
            % Compute fluctuating velocity components
            %
            % Args:
            %     obj (VelocityField): VelocityField instance with fields (u, v, w) containing the velocity components
            %
            % Optional Name-Pair Args:
            %     * density (float): Density field
            %     * weighted (bool): Flag to indicate if the velocity should be weighted by sqrt(density)
            %
            % Returns:
            %     velocity (VelocityField): Struct with fields (u, v, w) containing the fluctuating velocity components (fluctuations)
            %
            % Notes:
            %     * The velocity field is assumed to be periodic in all three dimensions.
            %     * By default, the velocity is weighted by density, if provided.
            %
            % Shortcuts:
            %     * velocity = getFluctuations(obj, rho);
            %
            % Examples:
            %     * velocity = getFluctuations(obj);
            %     * velocity = getFluctuations(obj, rho);
            %     * velocity = getFluctuations(obj, 'density', rho);
            %     * velocity = getFluctuations(obj, 'density', rho, 'weighted', true);
            
            % Default
            defaultDensity = [];
            defaultWeighted = true;

            % Shortcut: getFluctuations(obj, rho)
            if isscalar(varargin)
                varargin = {'density', varargin{1}};
            end

            % Parse input arguments
            p = inputParser;
            addParameter(p, 'density', defaultDensity, @(x) isnumeric(x));
            addParameter(p, 'weighted', defaultWeighted, @(x) islogical(x));
            parse(p, varargin{:});
            
            rho = p.Results.density;
            FLAG_WEIGHTED = p.Results.weighted;

            % Shallow copy of the velocity field
            velocity = obj.copy();

            % For compressible flows
            if ~isempty(rho)
                rhoMean = mean(rho, 'all');
                rhou = mean(rho .* velocity.u, 'all') / rhoMean;
                rhov = mean(rho .* velocity.v, 'all') / rhoMean;
                rhow = mean(rho .* velocity.w, 'all') / rhoMean;

                velocity.u = velocity.u - rhou;
                velocity.v = velocity.v - rhov;
                velocity.w = velocity.w - rhow;

                if FLAG_WEIGHTED
                    velocity.u = sqrt(rho) .* velocity.u;
                    velocity.v = sqrt(rho) .* velocity.v;
                    velocity.w = sqrt(rho) .* velocity.w;
                end

                return
            end

            % For incompressible flows
            velocity.u = velocity.u - mean(velocity.u, 'all');
            velocity.v = velocity.v - mean(velocity.v, 'all');
            velocity.w = velocity.w - mean(velocity.w, 'all');
        end

        function velocity = getCrossplaneFluctuations(obj, axisType, varargin)
            % Compute fluctuating velocity components
            %
            % Args:
            %     obj (VelocityField): VelocityField instance with fields (u, v, w) containing the velocity components
            %     axis (char): Axis for cross-plane averaging ('x', 'y', or 'z')
            %
            % Optional Name-Pair Args:
            %     * density (float): Density field
            %     * weighted (bool): Flag to indicate if the velocity should be weighted by sqrt(density)
            %
            % Returns:
            %     velocity (VelocityField): Struct with fields (u, v, w) containing the fluctuating velocity components (fluctuations)
            %
            % Notes:
            %     * The velocity field is assumed to be periodic in the direction perpendicular to the axis.
            %     * By default, the velocity is weighted by density, if provided.
            %
            % Examples:
            %     * velocity = getCrossplaneFluctuations(obj, 'x');
            %     * velocity = getCrossplaneFluctuations(obj, 'x', 'density', rho);
            %     * velocity = getCrossplaneFluctuations(obj, 'x', 'density', rho, 'weighted', true);

            % Default
            defaultDensity = [];
            defaultWeighted = true;

            % Parse input arguments
            p = inputParser;
            addRequired(p, 'axis', @(x) ischar(x) && ismember(lower(x), {'x', 'y', 'z'}))
            addParameter(p, 'density', defaultDensity, @(x) isnumeric(x));
            addParameter(p, 'weighted', defaultWeighted, @(x) islogical(x));
            parse(p, axisType, varargin{:});
            
            axisType = p.Results.axis;
            rho = p.Results.density;
            FLAG_WEIGHTED = p.Results.weighted;

            % Shallow copy of the velocity field
            velocity = obj.copy();

            % Map axis to slicing dimension
            sliceDim = find('xyz' == axisType);
            dims = 1:3;
            homogeneousDims = setdiff(dims, sliceDim);
            
            % For compressible flows
            if ~isempty(rho)
                rhoMean = mean(rho, homogeneousDims);
                rhou = mean(rho .* velocity.u, homogeneousDims) ./ rhoMean;
                rhov = mean(rho .* velocity.v, homogeneousDims) ./ rhoMean;
                rhow = mean(rho .* velocity.w, homogeneousDims) ./ rhoMean;

                velocity.u = velocity.u - rhou;
                velocity.v = velocity.v - rhov;
                velocity.w = velocity.w - rhow;

                if FLAG_WEIGHTED
                    velocity.u = sqrt(rho) .* velocity.u;
                    velocity.v = sqrt(rho) .* velocity.v;
                    velocity.w = sqrt(rho) .* velocity.w;
                end

                return
            end

            % For incompressible flows
            velocity.u = velocity.u - mean(velocity.u, homogeneousDims);
            velocity.v = velocity.v - mean(velocity.v, homogeneousDims);
            velocity.w = velocity.w - mean(velocity.w, homogeneousDims);
        end

        function velocityRMS = getVelocityRMS(obj)
            % Compute root mean square (rms) of the velocity field
            %
            % Args:
            %     obj (VelocityField): VelocityField object
            %
            % Returns:
            %     velocityRMS (float): Root mean square of the velocity field

            velocityRMS = sqrt( mean(obj.u.^2 + obj.v.^2 + obj.w.^2, 'all') / 3);
        end

        function Mt = getTurbulentMachNumberAverage(obj, soundVelocityAverage)
            % Compute mean turbulent Mach number
            %
            % Args:
            %     obj (VelocityField): VelocityField object
            %     soundVelocityAverage (float): Mean sound velocity
            %
            % Returns:
            %     Mt (float): Mean turbulent Mach number

            velocityRMS = obj.getVelocityRMS;
            Mt = sqrt(3) * velocityRMS / soundVelocityAverage;
        end

        function K = getTurbulentKineticEnergy(obj, varargin)
            % Compute Turbulent Kinetic Energy (TKE)
            %
            % Args:
            %     obj (VelocityField): VelocityField object
            %
            % Optional Args:
            %     * rho (float): 3D array with the density field
            %
            % Returns:
            %     K (float): 3D array with the TKE
            %
            % Example:
            %     K = getTurbulentKineticEnergy(u, v, w, rho);
            
            % Compressible flows
            if nargin > 1
                rho = varargin{1};
                K = 0.5 * mean(rho .* (obj.u.^2 + obj.v.^2 + obj.w.^2), 'all');
                return
            end

            % Incompressible flows
            K = 0.5 * mean((obj.u.^2 + obj.v.^2 + obj.w.^2), 'all');
        end

        function dissipation = getDissipation(obj, temperature, dynamicViscosity, varargin)
            % Compute dissipation energy
            %
            % Args:
            %     obj (VelocityField): VelocityField object
            %     temperature (float): 3D array with the temperature field
            %     dynamicViscosity (float): 3D array with the dynamic viscosity field
            %
            % Optional Args:
            %
            %
            % Returns:
            %     dissipation (float): 

            error('Method not implemented');
        end

        function omega_mag = getVorticity(obj)
            % Compute vorticity using FFT-based spectral differentiation
            %
            % Args:
            %     obj (VelocityField): VelocityField object with fields u, v, w
            %
            % Returns:
            %     omega_mag (float): 3D Array with the vorticity field
            %
            % Example:
            %     omega_mag = getVorticity(velocity)
        
            % Check if the grid is defined
            checkGridDefined(obj);

            % Get velocity FFTs
            [U, V, W] = getFFT(obj);
        
            % Get wave numbers
            [KX, KY, KZ] = getAngularWaveNumbers(obj);
        
            % Compute vorticity components
            omega_x = ifftn(1i * (KY .* W - KZ .* V), 'symmetric');
            omega_y = ifftn(1i * (KZ .* U - KX .* W), 'symmetric');
            omega_z = ifftn(1i * (KX .* V - KY .* U), 'symmetric');
        
            % Compute vorticity magnitude
            omega_mag = sqrt(omega_x.^2 + omega_y.^2 + omega_z.^2);
        end

        function div_u = getDivergence(obj)
            % Compute both divergence using FFT-based spectral differentiation
            %
            % Args:
            %     obj (VelocityField): VelocityField object with fields u, v, w
            %
            % Returns:
            %     div_u (float): 3D Array with the divergence of the velocity field
            %
            % Example:
            %     div_u = getDivergence(velocity)
        
            % Check if the grid is defined
            checkGridDefined(obj);

            % Get velocity FFTs
            [U, V, W] = getFFT(obj);
        
            % Get wave numbers
            [KX, KY, KZ] = getAngularWaveNumbers(obj);
        
            % Compute divergence
            div_u = ifftn(1i * (KX .* U + KY .* V + KZ .* W), 'symmetric');
        end

        function [omega_mag, div_u] = getVorticityDivergence(obj)
            % Compute both vorticity and divergence using FFT-based spectral differentiation
            %
            % Args:
            %     obj (VelocityField): VelocityField object with fields u, v, w
            %
            % Returns:
            %     Tuple containing:
            %
            %     * omega_mag (float): 3D Array with the vorticity field
            %     * div_u (float): 3D Array with the divergence of the velocity field
            %
            % Example:
            %     [omega_mag, div_u] = getVorticityDivergence(velocity)

            % Check if the grid is defined
            checkGridDefined(obj);

            % Get velocity FFTs
            [U, V, W] = getFFT(obj);
        
            % Get wave numbers
            [KX, KY, KZ] = getAngularWaveNumbers(obj);
        
            % Compute divergence
            div_u = ifftn(1i * (KX .* U + KY .* V + KZ .* W), 'symmetric');
        
            % Compute vorticity components
            omega_x = ifftn(1i * (KY .* W - KZ .* V), 'symmetric');
            omega_y = ifftn(1i * (KZ .* U - KX .* W), 'symmetric');
            omega_z = ifftn(1i * (KX .* V - KY .* U), 'symmetric');
        
            % Compute vorticity magnitude
            omega_mag = sqrt(omega_x.^2 + omega_y.^2 + omega_z.^2);
        end

        function magnitudeField = pointwiseMagnitude(obj)
            % Compute the pointwise magnitude of the velocity field
            %
            % Args:
            %     obj (VelocityField): VelocityField object
            %
            % Returns:
            %     magnitudeField (float): 3D matrix of velocity magnitudes
            %
            % Example:
            %     magnitudeField = pointwiseMagnitude(obj);

            magnitudeField = sqrt(obj.u.^2 + obj.v.^2 + obj.w.^2);
        end
        
        function globalMagnitude = globalMagnitude(obj)
            % Compute the global magnitude (Frobenius norm) of the velocity field
            %
            % Args:
            %     obj (VelocityField): VelocityField object
            %
            % Returns:
            %     globalMagnitude (float): Scalar magnitude of the velocity field
            %
            % Example:
            %     globalMagnitude = globalMagnitude(obj);

            globalMagnitude = sqrt(sum(obj.u(:).^2 + obj.v(:).^2 + obj.w(:).^2));
        end

        function obj = normalize(obj)
            % Normalize the velocity field by its magnitude
            %
            % Returns:
            %     obj: Normalized velocity field
            
            magnitude = obj.pointwiseMagnitude();
            obj.u = obj.u ./ magnitude;
            obj.v = obj.v ./ magnitude;
            obj.w = obj.w ./ magnitude;
        end

        function obj = scale(obj, factor)
            % Scale the velocity field by a factor
            %
            % Args:
            %     factor: Scaling factor
            %
            % Returns:
            %     obj: Scaled velocity field
            
            obj.u = obj.u * factor;
            obj.v = obj.v * factor;
            obj.w = obj.w * factor;
        end

        function [U, V, W] = getFFT(obj)
            % Get velocity FFTs
            %
            % Args:
            %     obj (VelocityField): VelocityField object with fields u, v, w
            %
            % Returns:
            %     Tuple containing:
            %
            %     * U (float): FFT of the u component
            %     * V (float): FFT of the v component
            %     * W (float): FFT of the w component
            %
            % Example:
            %     [U, V, W] = getFFT(obj)
            
            U = fftn(obj.u);
            V = fftn(obj.v);
            W = fftn(obj.w);
        end

        function ax = plotVorticity(obj, axisType, slices, varargin)
            % Plot vorticity magnitude
            %
            % Args:
            %     obj (VelocityField): VelocityField object
            %     axisType (char): Axis for slicing ('x', 'y', or 'z')
            %     slices (vector): Slices to plot
            % 
            % Optional Args:
            %     * omega_mag (float): 3D array with the vorticity magnitude
            %
            % Returns:
            %     ax: Axes handle of the plot
            
            % Import packages
            import combustiontoolbox.utils.extensions.brewermap

            % Check if omega_mag is provided
            if nargin < 4
                omega_mag = obj.getVorticity();
            else
                omega_mag = varargin{1};
            end

            % Clip at 99th percentile to avoid spurious peak values
            clip = prctile(omega_mag(:), 99);

            % Set color limits
            climits = [0, clip];
            
            % Plot 2D slices
            ax = obj.plotContour(omega_mag, axisType, slices, '$\nabla \times \mathbf{u}$', climits);
        end

        function ax = plotDivergence(obj, axisType, slices, varargin)
            % Plot divergence
            %
            % Args:
            %     obj (VelocityField): VelocityField object
            %     axisType (char): Axis for slicing ('x', 'y', or 'z')
            %     slices (vector): Slices to plot
            %
            % Optional Args:
            %     * div_u (float): 3D array with the divergence field
            %
            % Returns:
            %     ax: Axes handle of the plot
            
            % Check if div_u is provided
            if nargin < 4
                div_u = obj.getDivergence();
            else
                div_u = varargin{1};
            end

            % Normalize colormap
            absDiv = abs(div_u(:));
            clip = prctile(absDiv, 99);

            % Set color limits
            climits = [-clip, clip];
            
            % Plot 2D slices
            ax = obj.plotContour(div_u, axisType, slices, '$\nabla \cdot \mathbf{u}$', climits);
        end

    end
    
    methods (Access = private)

        function [Lx, Ly, Lz] = getDomainLengths(obj)
            % Get the lengths of the domain in each direction
            %
            % Args:
            %     obj (VelocityField): VelocityField object
            %
            % Returns:
            %     Tuple containing:
            %
            %     * Lx (float): Length in the x-direction
            %     * Ly (float): Length in the y-direction
            %     * Lz (float): Length in the z-direction
            
            % Check if the grid is defined
            checkGridDefined(obj);
            
            % Definitions
            Lx = obj.x(end) - obj.x(1);
            Ly = obj.y(end) - obj.y(1);
            Lz = obj.z(end) - obj.z(1);
        end

        function [KX, KY, KZ] = getAngularWaveNumbers(obj)
            % Compute angular wave numbers for FFT
            %
            % Args:
            %     sz (float): Size of the 3D array
            %
            % Returns:
            %     Tuple containing:
            %
            %     * KX (float): 3D array with the angular wave number in the x-direction
            %     * KY (float): 3D array with the angular wave number in the y-direction
            %     * KZ (float): 3D array with the angular wave number in the z-direction
            %
            % Example:
            %     [KX, KY, KZ] = getAngularWaveNumbers(velocity)

            % Check if the grid is defined
            checkGridDefined(obj);

            % Definitions
            sz = obj.size;
            [Lx, Ly, Lz] = getDomainLengths(obj);

            % Compute wave numbers
            kx = 2 * pi * ifftshift( -floor(sz(1)/2):ceil(sz(1)/2)-1 ) / Lx;
            ky = 2 * pi * ifftshift( -floor(sz(2)/2):ceil(sz(2)/2)-1 ) / Ly;
            kz = 2 * pi * ifftshift( -floor(sz(3)/2):ceil(sz(3)/2)-1 ) / Lz;

            % Create meshgrid for wave numbers
            [KX, KY, KZ] = ndgrid(kx, ky, kz);
        end
        
        function checkGridDefined(obj)
            % Check if the grid is defined
            %
            % Args:
            %     obj (VelocityField): VelocityField object
            %
            % Raises:
            %     error: If the grid is not defined
            
            if obj.gridChecked
                return
            end

            if isempty(obj.x) || isempty(obj.y) || isempty(obj.z)
                obj.x = linspace(0, obj.defaultLenght, obj.size(1));
                obj.y = linspace(0, obj.defaultLenght, obj.size(2));
                obj.z = linspace(0, obj.defaultLenght, obj.size(3));

                warning('VelocityField:GridDefaultUsed', ...
                        ['Grid coordinates (x, y, z) were not defined. ', ...
                         'Default uniform grids in [0, 2π] were used. ', ...
                         'It is recommended to provide physical grid coordinates.']);
            end

            % Update flag
            obj.gridChecked = true;
        end

    end

    methods (Access = private, Static)

        function ax = plotContour(field, axisType, slices, cLabel, climits)
            % Plot contour of a 2D slice of the velocity field
            %
            % Args:
            %     field (float): 3D array with the field to plot
            %     axisType (char): Axis for slicing ('x', 'y', or 'z')
            %     slices (vector): Slices to plot
            %     cLabel (string): Label for the colorbar
            %     climits (vector): Color limits for the plot
            %
            % Returns:
            %     ax: Axes handle of the plot
            %
            % Examples:
            %     * ax = VelocityField.plotContour(omega_mag, 'x', 1:10, '$\nabla \times \mathbf{u}$', [0, 10]);
            %     * ax = VelocityField.plotContour(div_u, 'x', 1:10, '$\nabla \cdot \mathbf{u}$, [-3, 3]);
            
             % Import packages
            import combustiontoolbox.utils.extensions.brewermap

            % Definitions
            numSlices = length(slices);
            map = brewermap(40, 'spectral');
            map = flip(map);
            fps = 1/300;

            % Map axis to slicing dimension
            sliceDim = find('xyz' == axisType);
            dims = 1:3;
            plotDims = setdiff(dims, sliceDim);

            % Permute dimensions to bring slice dimension last
            permutedField = permute(field, [sliceDim, plotDims]);
            
            ax = figure;

            for i = 1:numSlices
                slice = slices(i);
                contourf(squeeze(permutedField(slice, :, :)), 'LineColor', 'none');
                % imagesc(squeeze(permutedField(slice, :, :)));
                title( sprintf('Slice %d', slice), 'Interpreter', 'latex', 'FontSize', 20);
                axis equal off;
                clim(climits);
                colormap(map);

                % Add and customize colorbar
                cb = colorbar;
                cb.Label.Interpreter = 'latex';
                cb.Label.String = cLabel;
                cb.Label.FontSize = 18;
                cb.TickLabelInterpreter = 'latex';
                cb.FontSize = 16;

                pause(fps);
            end

        end

    end

end