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

        function obj = plus(obj, obj2)
            % Overload the plus operator to add the velocity field from two VelocityField objects
            %
            % Args:
            %     obj (VelocityField): VelocityField object
            %     obj2 (VelocityField): VelocityField object
            %
            % Returns:
            %     obj (VelocityField): VelocityField object 

            obj.u = obj.u + obj2.u;
            obj.v = obj.v + obj2.v;
            obj.w = obj.w + obj2.w;            
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
            
            magnitude = obj.computeMagnitude();
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

    end

end