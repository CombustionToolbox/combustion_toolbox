classdef VelocityField < handle
    % This :mat:func:`VelocityField` class provides methods for computing properties of a velocity field.
    
    properties
        u % x-component of the velocity field
        v % y-component of the velocity field
        w % z-component of the velocity field
    end

    methods
        
        function obj = VelocityField(u, v, w)
            % Constructor for VelocityField
            %
            % Args:
            %     u (float): x-component of velocity
            %     v (float): y-component of velocity
            %     w (float): z-component of velocity
            %
            % Returns:
            %     obj (VelocityField): VelocityField object
            
            % Parse inputs
            p = inputParser;
            addRequired(p, 'u', @isnumeric);
            addRequired(p, 'v', @isnumeric);
            addRequired(p, 'w', @isnumeric);
            parse(p, u, v, w);

            % Assign properties
            obj.u = p.Results.u;
            obj.v = p.Results.v;
            obj.w = p.Results.w;
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