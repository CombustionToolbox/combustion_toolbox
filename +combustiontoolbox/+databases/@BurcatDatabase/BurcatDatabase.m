classdef BurcatDatabase < combustiontoolbox.databases.Database
    
    methods (Access = public)
        
        function obj = BurcatDatabase(varargin)
            % Constructor 

            % Call superclass constructor
            obj@combustiontoolbox.databases.Database('name', 'Burcat', 'temperatureReference', 298.15, varargin{:});
        end
    
    end

end