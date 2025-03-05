classdef MixtureConfig < handle
    % The :mat:class:`MixtureConfig` class is used to store configuration settings 
    % for a :mat:class:`Mixture` object.
    %
    % The :mat:class:`MixtureConfig` object can be initialized as follows: ::
    %
    %      config = MixtureConfig('mintolDisplay', 1e-6, ...
    %                             'compositionUnits', 'molar fraction', ...
    %                             'FLAG_COMPACT', true)
    %
    % This creates an instance of the `MixtureConfig` class with the specified configuration.
    %
    % See also: :mat:class:`Mixture`
    
    properties
        mintolDisplay = 1e-6                % Minimum tolerance to display species
        compositionUnits = 'molar fraction' % Options: 'mol', 'molar fraction', or 'mass fraction'
        FLAG_COMPACT = true                 % FLAG to print composition in a compact format
    end
    
    methods
        function obj = MixtureConfig(varargin)
            % MixtureConfig contructor
            %
            % This function accepts name-value pair arguments to override default property values.
            %
            % Examples:
            %     * config = MixtureConfig();
            %     * config = MixtureConfig('mintolDisplay', 1e-6);
            %     * config = MixtureConfig('compositionUnits', 'molar fraction');
            %     * config = MixtureConfig('FLAG_COMPACT', true);
            %     * config = MixtureConfig('mintolDisplay', 1e-6, 'compositionUnits', 'molar fraction', 'FLAG_COMPACT', true);
            
            % Definitions
            validUnits = {'mol', 'molar fraction', 'mass fraction'};

            % Parse inputs
            ip = inputParser;
            % Validate mintolDisplay: must be a positive scalar numeric value
            addParameter(ip, 'mintolDisplay', obj.mintolDisplay, @(x) isnumeric(x) && isscalar(x) && x > 0);
            addParameter(ip, 'compositionUnits', obj.compositionUnits, @(x) (isstring(x) || ischar(x)) && any(strcmp(x, validUnits)));
            addParameter(ip, 'FLAG_COMPACT', obj.FLAG_COMPACT, @(x) islogical(x) && isscalar(x));
            parse(ip, varargin{:});
            
            % Set properties
            obj.mintolDisplay = ip.Results.mintolDisplay;
            obj.compositionUnits = ip.Results.compositionUnits;
            obj.FLAG_COMPACT = ip.Results.FLAG_COMPACT;
        end
    end
end