classdef Integrator < handle
    % Integrator - Abstract base class for numerical integration
    %
    % This class provides a framework for numerical integration methods.
    % It is intended to be subclassed for specific integration techniques.
    %
    % Example:
    %   integrator = Integrator();
    %   result = integrator.integrate(@(x) x.^2, 0, 1);
    
    properties
        tolRelative  = 1e-6     % Relative tolerance for integrals
        tolAbsolute  = 1e-10    % Absolute tolerance for integrals
    end
    
    methods
        
        function obj = Integrator(varargin)
            % Constructor

            % Parse input arguments
            p = inputParser;
            addParameter(p, 'tolRelative', obj.tolRelative, @(x) isnumeric(x) && isscalar(x));
            addParameter(p, 'tolAbsolute', obj.tolAbsolute, @(x) isnumeric(x) && isscalar(x));
            parse(p, varargin{:});
            
            % Set properties
            obj.tolRelative = p.Results.tolRelative;
            obj.tolAbsolute = p.Results.tolAbsolute;
        end
        
        function value = integrate(obj, fun, a, b, varargin)
            % Integrate a function fun from a to b using numerical integration
            %
            % Args:
            %     fun (function): Function handle to be integrated
            %     a (float): Lower limit of integration
            %     b (float): Upper limit of integration
            %
            % Returns:
            %     value(float): result of the integration
            
            value = integral(fun, a, b, 'RelTol', obj.tolRelative, 'AbsTol', obj.tolAbsolute, varargin{:});
        end
        
    end
    
end