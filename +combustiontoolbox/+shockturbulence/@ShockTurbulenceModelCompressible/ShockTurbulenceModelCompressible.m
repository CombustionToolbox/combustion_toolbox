classdef ShockTurbulenceModelCompressible < combustiontoolbox.shockturbulence.ShockTurbulenceModel
    % The :mat:class:`ShockTurbulenceModelCompressible` class characterizes 
    % turbulence amplification across a shock wave interacting with weak turbulence
    % comprised of vortical-entropic and acoustic disturbances using linear theory.
    %
    % These models are based on our previous theoretical work [1]
    % and have been extended to multi-component mixtures [2, 3] using the
    % Combustion Toolbox [4].
    %
    % References:
    %     [1] Huete, C., Cuadra, A., Vera, M., Urzay, & J. (2021). Thermochemical
    %         effects on hypersonic shock waves interacting with weak turbulence.
    %         Physics of Fluids 33, 086111 (featured article). DOI: 10.1063/5.0059948.
    %
    %     [2] Cuadra, A., Vera, M., Di Renzo, M., & Huete, C. (2023). Linear Theory
    %         of Hypersonic Shocks Interacting with Turbulence in Air. In 2023 AIAA
    %         SciTech Forum, National Harbor, USA. DOI: 10.2514/6.2023-0075.
    %
    %     [3] Cuadra, A., Di Renzo, M., Hoste, J. J. O., Williams, C. T., Vera, M., & Huete, C. (2025).
    %         Review of shock-turbulence interaction with a focus on hypersonic flow. Physics of Fluids, 37(4).
    %         DOI: 10.1063/5.0255816.
    %
    %     [4] Cuadra, A., Huete, C., Vera, M. (2022). Combustion Toolbox:
    %         A MATLAB-GUI based open-source tool for solving gaseous
    %         combustion problems. Zenodo. DOI: 10.5281/zenodo.5554911.
    
    properties
        shockTurbulenceModelVorticalEntropic
        shockTurbulenceModelAcoustic
    end

    methods

        function obj = ShockTurbulenceModelCompressible(varargin)
            % Constructor
            %
            % Optional Args:
            %   varargin (optional): key-value pairs to initialize the database
            %
            % Returns:
            %     obj (ShockTurbulenceModelCompressible): ShockTurbulenceModelCompressible object
            %
            % Examples:
            %     * shockTurbulenceModel = combustiontoolbox.shockturbulence.ShockTurbulenceModelCompressible();
            %     * shockTurbulenceModel = combustiontoolbox.shockturbulence.ShockTurbulenceModelCompressible('dimensions', 3);
            %     * shockTurbulenceModel = combustiontoolbox.shockturbulence.ShockTurbulenceModelCompressible('tolIntegralRelative', 1e-6);
            %     * shockTurbulenceModel = combustiontoolbox.shockturbulence.ShockTurbulenceModelCompressible('tolIntegralAbsolute', 1e-10);

            % Set shock turbulence models
            obj.shockTurbulenceModelVorticalEntropic = combustiontoolbox.shockturbulence.ShockTurbulenceModelVorticalEntropic(varargin{:});
            obj.shockTurbulenceModelAcoustic = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic(varargin{:});
        end
            
        function averages = getAverages(obj, R, M2, Gammas, Gammas1, Gammas3, beta, eta, chi)
            % getAverages Compute the amplification ratios.
            %
            %   results = getAverages(R, M2, Gammas, Gammas1, Gammas3, beta, eta, chi)
            %
            % Args:
            %     obj (ShockTurbulenceModelCompressible): ShockTurbulenceModelCompressible object
            %     R (float): Density ratio rho_2 / rho_1
            %     M2 (float): Post-shock Mach number
            %     Gammas (float): Inverse normalized Hugonito slope, see Eq. (22) in [1] (Eq. (4) in [2])
            %     Gammas1 (float): Inverse normalized Hugonito slope
            %     Gammas3 (float): Inverse normalized Hugonito slope
            %     beta (float): Ratio of speed of sound in the post-shock state to the pre-shock state
            %     eta (float): square ratio of vortical to acoustic disturbances
            %     chi (float): mean correlation coefficient between the vortical and entropic modes
            %
            % Returns:
            %     averages (struct): Structure with computed averages (K, L, T, etc.)
            %
            % Example:
            %     averages = getAverages(ShockTurbulenceModelCompressible(), R, M2, Gammas, Gammas1, Gammas3, beta, eta, chi);
            
            % Parse input arguments
            p = inputParser;
            addRequired(p, 'R', @(x) isnumeric(x));
            addRequired(p, 'M2', @(x) isnumeric(x));
            addRequired(p, 'Gammas', @(x) isnumeric(x));
            addRequired(p, 'Gammas1', @(x) isnumeric(x));
            addRequired(p, 'Gammas3', @(x) isnumeric(x));
            addRequired(p, 'beta', @(x) isnumeric(x));
            addRequired(p, 'eta', @(x) isnumeric(x));
            addRequired(p, 'chi', @(x) isnumeric(x));
            parse(p, R, M2, Gammas, Gammas1, Gammas3, beta, eta, chi);

            % Set properties
            R = p.Results.R;
            M2 = p.Results.M2;
            Gammas = p.Results.Gammas;
            Gammas1 = p.Results.Gammas1;
            Gammas3 = p.Results.Gammas3;
            beta = p.Results.beta;
            eta = p.Results.eta;
            chi = p.Results.chi;

            % Compute results assuming vortical-entropic fluctuations
            averagesVorticalEntropic = getAverages(obj.shockTurbulenceModelVorticalEntropic, R, M2, Gammas, Gammas1, beta, chi);

            % Compute results assuming acoustic fluctuations
            averagesAcoustic = getAverages(obj.shockTurbulenceModelAcoustic, R, M2, Gammas, Gammas1, Gammas3, beta);

            % Compute post-shock turbulence statistics by supersposition of the vortical-entropic and acoustic solutions
            averages.K = (averagesVorticalEntropic.K + eta .* averagesAcoustic.K) ./ (1 + eta);
            averages.R11 = (averagesVorticalEntropic.R11 + eta .* averagesAcoustic.R11) ./ (1 + eta);
            averages.RTT = (averagesVorticalEntropic.RTT + eta .* averagesAcoustic.RTT) ./ (1 + eta);
            averages.Ka = (averagesVorticalEntropic.Ka + eta .* averagesAcoustic.Ka) ./ (1 + eta);
            averages.Kr = (averagesVorticalEntropic.Kr + eta .* averagesAcoustic.Kr) ./ (1 + eta);
            averages.enstrophy = (averagesVorticalEntropic.enstrophy + eta .* averagesAcoustic.enstrophy) ./ (1 + eta);
            averages.enstrophyTT = (averagesVorticalEntropic.enstrophyTT + eta .* averagesAcoustic.enstrophyTT) ./ (1 + eta);
        end

        function obj = setPDF(obj)
            % Define the probability density function (PDF) used to weigh different
            % wavenumber contributions in the post-shock turbulence model
            %
            % Args:
            %     obj (ShockTurbulenceModelCompressible): ShockTurbulenceModelCompressible object
            %
            % Returns:
            %     obj (ShockTurbulenceModelCompressible): Updated ShockTurbulenceModelCompressible object with PDF set
            
            warning('setPDF method not needed.');
        end

    end

end