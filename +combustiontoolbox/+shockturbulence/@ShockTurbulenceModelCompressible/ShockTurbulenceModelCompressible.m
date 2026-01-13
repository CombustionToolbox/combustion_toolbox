classdef ShockTurbulenceModelCompressible < combustiontoolbox.shockturbulence.ShockTurbulenceModel
    % The :mat:func:`ShockTurbulenceModelCompressible` class characterizes 
    % turbulence amplification across a shock wave interacting with weak turbulence
    % comprised of vortical-entropic and acoustic disturbances using linear theory.
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
            
        function [averages, mixArray1, mixArray2] = getAverages(obj, jumpConditions, mixArray1, mixArray2)
            % Compute the post-shock turbulence statistics by superposing the contributions of acoustic and vortical-entropic fluctuations
            %
            % Args:
            %     obj (ShockTurbulenceModelCompressible): ShockTurbulenceModelCompressible object
            %     jumpConditions (struct): Structure with jump conditions across the shock wave
            %     mixArray1 (Mixture): Pre-shock Mixture objects
            %     mixArray2 (Mixture): Post-shock Mixture objects
            %
            % Returns:
            %     Tuple containing:
            %
            %     * averages (struct): Structure with averages of post-shock turbulence statistics (e.g., Reynolds stresses, turbulent kinetic energy, enstrophy, etc.)
            %     * mixArray1 (Mixture): Pre-shock Mixture objects
            %     * mixArray2 (Mixture): Post-shock Mixture objects
            %
            % Example:
            %     [averages, mixArray1, mixArray2] = getAverages(ShockTurbulenceModelCompressible(), jumpConditions, mixArray1, mixArray2);

            % Definitions
            beta = jumpConditions.beta;
            eta = [mixArray1.eta];
            etaVorticity = [mixArray1.etaVorticity];

            % Compute results assuming vortical-entropic fluctuations
            averagesVorticalEntropic = getAverages(obj.shockTurbulenceModelVorticalEntropic, jumpConditions, mixArray1, mixArray2);

            % Compute results assuming acoustic fluctuations
            averagesAcoustic = getAverages(obj.shockTurbulenceModelAcoustic, jumpConditions, mixArray1, mixArray2);

            % Compute post-shock turbulence statistics by supersposition of the vortical-entropic and acoustic solutions
            averages.K = (averagesVorticalEntropic.K + eta .* averagesAcoustic.K) ./ (1 + eta);
            averages.R11 = (averagesVorticalEntropic.R11 + eta .* averagesAcoustic.R11) ./ (1 + eta);
            averages.RTT = (averagesVorticalEntropic.RTT + eta .* averagesAcoustic.RTT) ./ (1 + eta);
            averages.Ka = (averagesVorticalEntropic.Ka + eta .* averagesAcoustic.Ka) ./ (1 + eta);
            averages.Kr = (averagesVorticalEntropic.Kr + eta .* averagesAcoustic.Kr) ./ (1 + eta);
            averages.R11a = (averagesVorticalEntropic.R11a + eta .* averagesAcoustic.R11a) ./ (1 + eta);
            averages.R11r = (averagesVorticalEntropic.R11r + eta .* averagesAcoustic.R11r) ./ (1 + eta);
            averages.RTTa = (averagesVorticalEntropic.RTTa + eta .* averagesAcoustic.RTTa) ./ (1 + eta);
            averages.RTTr = (averagesVorticalEntropic.RTTr + eta .* averagesAcoustic.RTTr) ./ (1 + eta);
            averages.enstrophy = (averagesVorticalEntropic.enstrophy + etaVorticity .* beta.^2 .* averagesAcoustic.enstrophy);
            averages.enstrophyTT = (averagesVorticalEntropic.enstrophyTT + etaVorticity .* beta.^2 .* averagesAcoustic.enstrophyTT);

            % Get Kolmogorov length scale ratio across the shock
            averages.kolmogorovLengthRatio = obj.getKolmogorovLength(averages, mixArray1, mixArray2);
            
            % Set the post-shock turbulence statistics in the Mixture array
            mixArray2 = obj.setAverages2MixArray(averages, mixArray2);
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