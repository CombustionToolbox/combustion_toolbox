classdef (Abstract) ShockTurbulenceModel < handle
    % The :mat:func:`ShockTurbulenceModel` abstract class defines the interface 
    % for modeling turbulence amplification across a shock wave interacting with 
    % weak turbulence using linear interaction analysis. It serves as the base class 
    % for specific shock–turbulence interaction models.
    %
    % Uptream turbulence can be comprised of the following type of disturbances:
    % 	* vortical (shockTurbulenceModelVortical)
    %   * acoustic (shockTurbulenceModelAcoustic)
    %   * vortical + entropic (shockTurbulenceModelVorticalEntropic)
    %   * vortical + entropic + acoustic (shockTurbulenceModelCompressible)
    % 
    % These models are based on our previous theoretical works [1-3]
    % and have been extended to multi-component mixtures [4-7] using the
    % Combustion Toolbox [8-9].
    %
    % References:
    %     [1] Huete, C., Cuadra, A., Vera, M., Urzay, & J. (2021). Thermochemical
    %         effects on hypersonic shock waves interacting with weak turbulence.
    %         Physics of Fluids 33, 086111 (featured article). DOI: 10.1063/5.0059948.
    %
    %     [2] Huete, C., Velikovich, A. L., & Wouchuk, J. G. (2011). Analytical linear theory
    %         for the interaction of a planar shock wave with a two-or three-dimensional
    %         random isotropic density field. Physical Review E—Statistical, Nonlinear, and
    %         Soft Matter Physics, 83(5), 056320. DOI: 10.1103/PhysRevE.83.056320.
    %
    %     [3] Huete, C., Wouchuk, J. G., & Velikovich, A. L. (2012). Analytical linear theory
    %         for the interaction of a planar shock wave with a two-or three-dimensional random
    %         isotropic acoustic wave field. Physical Review E—Statistical, Nonlinear, and Soft
    %         Matter Physics, 85(2), 026312. DOI: 10.1063/5.0059948.
    %
    %     [4] Cuadra, A., Vera, M., Di Renzo, M., & Huete, C. (2023). Linear Theory
    %         of Hypersonic Shocks Interacting with Turbulence in Air. In 2023 AIAA
    %         SciTech Forum, National Harbor, USA. DOI: 10.2514/6.2023-0075.
    %
    %     [5] Cuadra, A., Williams, C. T., Di Renzo, M. & Huete, C. (2024). Compressibility
    %         and vibrational-excitation effects in hypersonic shock-turbulence interaction.
    %         Tech. Rep. Summer Program Proceedings, Center for Turbulence Research,
    %         Stanford University.
    %
    %     [6] Cuadra, A., Williams, C. T., Di Renzo, M., & Huete, C. The role of compressibility
    %         and vibrational-excitation in hypersonic shock–turbulence interactions.
    %         Journal of Fluid Mechanics (under review).
    %
    %     [7] Cuadra, A., Di Renzo, M., Hoste, J. J. O., Williams, C. T., Vera, M., & Huete, C. (2025).
    %         Review of shock-turbulence interaction with a focus on hypersonic flow. Physics of Fluids, 37(4).
    %         DOI: 10.1063/5.0255816.
    %
    %     [8] Cuadra, A., Huete, C., & Vera, M. (2026). Combustion Toolbox: An open-source
    %         thermochemical code for gas-and condensed-phase problems involving chemical equilibrium. 
    %         Computer Physics Communications 320, 110004. DOI:10.1016/j.cpc.2025.110004.
    %
    %     [9] Cuadra, A., Huete, C., Vera, M. (2022). Combustion Toolbox:
    %         A MATLAB-GUI based open-source tool for solving gaseous
    %         combustion problems. Zenodo. DOI: 10.5281/zenodo.5554911.

    properties
        problemType    % Type of problem ('compressible', 'vortical', 'vortical_entropic', 'acoustic')  
        averages       % A structure with computed averages
        dimensions     % Number of dimensions (2D or 3D turbulence)
        integrator     % Integrator object for numerical integration
        viscosityModel % Viscosity model ('powerlaw' or 'sutherland'). This is a temporal function and will be overriden with a specific TransportProperties class in future releases.
    end

    properties (Access = protected)
        pdf            % Probability density function (PDF) for the turbulence model
    end
    
    methods

        function obj = ShockTurbulenceModel(varargin)
            % Constructor for ShockTurbulenceModel
            
            % Default values
            defaultProblemType = 'vortical';
            defaultDimensions = 3;
            defaultIntegrator = combustiontoolbox.utils.math.Integrator();
            defaultViscosityModel = 'sutherland';

            % Parse inputs
            p = inputParser;
            addOptional(p, 'problemType', defaultProblemType, @(x) ischar(x) && any(strcmpi(x, {'compressible', 'vortical', 'vortical_entropic', 'acoustic'})));
            addParameter(p, 'dimensions', defaultDimensions, @(x) isnumeric(x) && isscalar(x) && (x >= 1 || x <= 3));
            addParameter(p, 'integrator', defaultIntegrator, @(x) isa(x, 'combustiontoolbox.utils.math.Integrator'));
            addParameter(p, 'viscosityModel', defaultViscosityModel, @(x) ischar(x) && any(strcmpi(x, {'powerlaw', 'sutherland'})));
            parse(p, varargin{:});
            
            % Set properties
            obj.problemType = p.Results.problemType;
            obj.dimensions = p.Results.dimensions;
            obj.integrator = p.Results.integrator;
            obj.viscosityModel = p.Results.viscosityModel;
        end

    end


    methods (Abstract)
        averages = getAverages(obj, varargin)
        obj = setPDF(obj)
    end
    
    methods (Access = protected)

        function value = integrate(obj, fun, a, b, varargin)
            % Integrate a function fun from a to b using the Integrator class
            %
            % Args:
            %     fun (function): Function handle to be integrated
            %     a (float): Lower limit of integration
            %     b (float): Upper limit of integration
            %
            % Returns:
            %     value(float): result of the integration
            
            value = obj.integrator.integrate(fun, a, b, varargin{:});
        end

    end

end
