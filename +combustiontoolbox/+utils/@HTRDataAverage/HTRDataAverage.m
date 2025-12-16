classdef HTRDataAverage < handle & combustiontoolbox.utils.HTRDataReader
    % The :mat:func:`HTRDataAverage` class provides methods to process averaged
    % data from the HTR solver [1] output files, which are stored in the
    % Hierarchical Data Format (HDF) file format.
    %
    % The Hypersonic Task-based Research (HTR) solver [1] is an open-source
    % exascale-oriented task-based multi-GPU high-order code for hypersonics
    % aerothermodynamics. The HTR solver is available at the GitHub repository [2].
    %
    % Note: The root-mean-square (rms) values from the HTR solver output files
    % are not the true rms values, they are the square of the rms values.
    %
    % References:
    %   [1] Di Renzo, M., Fu, L., and Urzay, J., HTR solver: An open-source
    %       exascale-oriented task-based multi-GPU high-order code for
    %       hypersonics aerothermodynamics, Comput. Phys. Commun, Vol. 255,
    %       2020, p. 107262
    %   [2] https://github.com/stanfordhpccenter/HTR-solver
    
    properties
        averages              % Struct holding computed averages and derived data
        filepath              % Path to the .hdf file
        axis = 'x'            % Axis for cross-plane averaging ('x', 'y', or 'z')
        k0                    % Energy-peak wavenumber
        eta                   % Ratio of dilatational to solenoidal turbulent kinetic energy
        tolShock = 1e-4       % Relative tolerance for shock detection
        indexPreShock         % Index of the pre-shock location
        indexShock            % Index of the mean shock location
        indexPostShock        % Index of the post-shock location
        indexFarfield         % Index of the far-field location to compare with LIA results
        plotConfig            % PlotConfig object
        plotYaxisScale = [0.9, 1.3] % Factor to scale the y-axis in the plot
        FLAG_SHIFT = true     % Flag to shift the grid coordinate to the shock location
    end

    properties (Access = private)
        rawData               % Struct holding raw data from the .hdf file
        weight                % Data weighting factor
        sliceDim              % Slicing dimension
        numProperties         % Number of properties in the raw data
        numProperties_rms     % Number of properties in the raw data with rms
        numProperties_favg    % Number of properties in the raw data with Favre average
        numProperties_frms    % Number of properties in the raw data with Favre rms
        listProperties        % List of properties in the raw data
        listProperties_rms    % List of properties in the raw data with rms
        listProperties_favg   % List of properties in the raw data with Favre average
        listProperties_frms   % List of properties in the raw data with Favre rms
        FLAG_K0 = false       % Flag to normalize the grid coordinate
    end

    methods

        function obj = HTRDataAverage(filepath, varargin)
            % Constructor for HTRDataAverage
            %
            % Args:
            %     filepath (char): Path to the .hdf file
            
            % Default
            defaultPlotConfig = combustiontoolbox.utils.display.PlotConfig();
            defaultPlotConfig.innerposition = [0.15 0.15 0.35 0.5];
            defaultPlotConfig.outerposition = [0.15 0.15 0.35 0.5];

            % Parse input arguments
            p = inputParser;
            addRequired(p, 'filepath', @ischar);
            addParameter(p, 'axis', obj.axis, @(x) ischar(x) && ismember(lower(x), {'x', 'y', 'z'}));
            addParameter(p, 'k0', [], @(x) isnumeric(x) && x > 0);
            addParameter(p, 'eta', [], @(x) isnumeric(x) && x >= 0);
            addParameter(p, 'tolShock', obj.tolShock, @(x) isnumeric(x) && x > 0);
            addParameter(p, 'plotConfig', defaultPlotConfig, @(x) isa(x, 'combustiontoolbox.utils.display.PlotConfig'));
            addParameter(p, 'FLAG_SHIFT', obj.FLAG_SHIFT, @islogical);
            parse(p, filepath, varargin{:});

            % Set properties
            obj.filepath = p.Results.filepath;
            obj.axis = p.Results.axis;
            obj.k0 = p.Results.k0;
            obj.eta = p.Results.eta;
            obj.plotConfig = p.Results.plotConfig;
            obj.tolShock = p.Results.tolShock;
            obj.FLAG_SHIFT = p.Results.FLAG_SHIFT;
            
            % Get listProperties
            obj.listProperties = {h5info(obj.filepath).Datasets.Name};
            obj.listProperties_rms = obj.listProperties(contains(obj.listProperties, '_rms'));
            obj.listProperties_favg = obj.listProperties(contains(obj.listProperties, '_favg'));
            obj.listProperties_frms = obj.listProperties(contains(obj.listProperties, '_frms'));

            % Definitions
            obj.numProperties = length(obj.listProperties);
            obj.numProperties_rms = length(obj.listProperties_rms);
            obj.numProperties_favg = length(obj.listProperties_favg);
            obj.numProperties_frms = length(obj.listProperties_frms);

            % Map axis to slicing dimension
            obj.sliceDim = find('xyz' == obj.axis);

            % Load raw data and apply weight
            load(obj);

            % Compute averages
            getAverages(obj);

            % Get grid coordinate in the axis for cross-plane averaging
            getGrid(obj);

            % Find shock location
            getShockLocation(obj);

            % Compute local turbulent Mach number
            getLocalTurbulentMachNumber(obj);

            % Compute dissipation rate
            getDissipation(obj);

            % Compute Kolmogorov length scale
            getKolmogorovLength(obj);

            % Compute turbulent kinetic energy (TKE) amplification ratio 
            getTurbulentKineticEnergy(obj);

            % Get far-field location
            getFarfieldLocation(obj);

            % Compute variances
            getVariances(obj);

            % Compute correlations
            getCorrelations(obj);

            % Print summary
            print(obj);
        end

        function load(obj, varargin)
            % Load raw data from the HTR solver output file

            % Check if filepath is provided
            if nargin > 1
                obj.filepath = varargin{1};
            end

            % Load data
            obj.weight(1, :) = obj.read(obj.filepath, 'weight');

            % Load raw data and apply weight
            for j = obj.numProperties:-1:1
                % Get property
                property = obj.listProperties{j};

                % Read data
                obj.rawData.(property) = obj.read(obj.filepath, property);

                % If the first dimension is not the direction, reshape
                sz = size(obj.rawData.(property));
                if sz(1) ~= min(sz)
                    obj.rawData.(property) = obj.rawData.(property)';
                end

                % Apply weight
                applyWeight(obj, property);
            end

        end

        function getAverages(obj)
            % Complet averaging of the raw data

            % Compute root mean square (rms)
            for property = obj.listProperties_rms
                property = property{:};
                obj.averages.(property) = getRMS(obj, property);
            end

            % Compute Favre root mean square (rms)
            for property = obj.listProperties_frms
                property = property{:};
                obj.averages.(property) = getFavreRMS(obj, property);
            end

            % Compute Favre velocity fluctuations
            obj.averages.velocity_frey(1, :) = obj.averages.velocity_frey(1, :) - obj.averages.velocity_favg(1, :) .* obj.averages.velocity_favg(2, :) ./ obj.averages.rho_avg;
            obj.averages.velocity_frey(2, :) = obj.averages.velocity_frey(2, :) - obj.averages.velocity_favg(1, :) .* obj.averages.velocity_favg(3, :) ./ obj.averages.rho_avg;
            obj.averages.velocity_frey(3, :) = obj.averages.velocity_frey(3, :) - obj.averages.velocity_favg(2, :) .* obj.averages.velocity_favg(3, :) ./ obj.averages.rho_avg;
        end

        function [x, dx] = getGrid(obj)
            % Get grid coordinate in the axis for cross-plane averaging and find shock location
            %
            % Args:
            %     obj (HTRDataAverage): HTRDataAverage object
            %
            % Returns:
            %     x (float): 1D array with the grid coordinate in the axis for cross-plane averaging
            %     dx (float): 1D array with the grid spacing in the axis for cross-plane averaging

            % Get grid coordinate in the axis for cross-plane averaging
            x = obj.averages.centerCoordinates(obj.sliceDim, :);

            % Normalize grid coordinate if k0 is provided
            if ~isempty(obj.k0)
                x = x * obj.k0;
                obj.FLAG_K0 = true;
            end

            % Compute grid spacing in the axis for cross-plane averaging
            dx = diff(x);

            % Assign properties
            obj.averages.x = x;
            obj.averages.dx = dx;
        end

        function [xShock, indexShock] = getShockLocation(obj)
            % Identify the mean shock location based on velocity gradients
            %
            % Args:
            %     obj (HTRDataAverage): HTRDataAverage object
            %
            % Returns:
            %     xShock (float): Shock location
            %     indexShock (float): Index of the mean shock location

            % Compute gradient of the velocity field
            dudx = diff(obj.averages.velocity_avg(obj.sliceDim, :)) ./ obj.averages.dx;

            % Find pre-shock shock index
            [~, indexPreShock] = find(abs(dudx) > obj.tolShock, 1);

            % Find shock index
            [~, indexShock] = max(-dudx);
            
            % Find post-shock index
            indexPostShock = find(dudx(indexShock:end) / max(dudx(indexShock:end)) > obj.tolShock, 1) + indexShock;
            
            % Compute shock location
            xShock = obj.averages.x(indexPostShock);

            % Shift grid coordinate to the shock location
            if obj.FLAG_SHIFT
                obj.averages.x = obj.averages.x - obj.averages.x(indexShock);
            end

            % Assign properties
            obj.indexPreShock = indexPreShock;
            obj.indexShock = indexShock;
            obj.indexPostShock = indexPostShock;
        end

        function Mt = getLocalTurbulentMachNumber(obj)
            % Compute local turbulent Mach number from the data
            %
            % Args:
            %     obj (HTRDataAverage): HTRDataAverage object
            %
            % Returns:
            %     Mt (float): Local turbulent Mach number

            % Compute local turbulent Mach number
            Mt = sqrt(sum(obj.averages.velocity_rms)) ./ obj.averages.SoS_avg;

            % Assign properties
            obj.averages.Mt = Mt;
        end

        function [dissipation, dissipationVolumetric] = getDissipation(obj)
            % Compute dissipation rate from the data
            %
            % Args:
            %     obj (HTRDataAverage): HTRDataAverage object
            %
            % Returns:
            %     dissipation (float): Dissipation rate

            % Compute volumetric dissipation rate
            dissipationVolumetric = obj.averages.tauGradU(1, :) + obj.averages.tauGradU(2, :) + obj.averages.tauGradU(3, :);

            % Compute dissipation rate per unit mass
            dissipation = dissipationVolumetric ./ obj.averages.rho_avg;

            % Assign properties
            obj.averages.dissipation = dissipation;
            obj.averages.dissipationVolumetric = dissipationVolumetric;
        end

        function [lengthKolmogorov, ratioGridLengthKolmogorov] = getKolmogorovLength(obj)
            % Compute Kolmogorov length scale from the data
            %
            % Args:
            %     obj (HTRDataAverage): HTRDataAverage object
            %
            % Returns:
            %     lengthKolmogorov (float): Kolmogorov length scale
            %     ratioGridLengthKolmogorov (float): Ratio of the grid length to the Kolmogorov length scale

            % Definitions
            nu = obj.averages.mu_avg ./ obj.averages.rho_avg;
            dissipation = abs( getDissipation(obj) );
            dx = obj.averages.dx;

            % Remove normalization if k0 is provided
            if obj.FLAG_K0
                dx  = dx / obj.k0;
            end

            % Compute Kolmogorov length scale
            lengthKolmogorov = (nu.^3 ./ dissipation).^(1/4);
            
            % Compute ratio of the grid length to the Kolmogorov length scale
            ratioGridLengthKolmogorov = dx ./ lengthKolmogorov(1:end-1);

            % Assign properties
            obj.averages.lengthKolmogorov = lengthKolmogorov;
            obj.averages.lengthKolmogorovk0 = lengthKolmogorov * obj.k0;
            obj.averages.ratioLengthKolmogorov = lengthKolmogorov / lengthKolmogorov(obj.indexPreShock);
            obj.averages.ratioGridLengthKolmogorov = ratioGridLengthKolmogorov;
            obj.averages.ratioKolmogorovGridLength = 1 ./ ratioGridLengthKolmogorov;
        end

        function [TKE_f, R11_f, R22_f, R33_f, RTT_f] = getTurbulentKineticEnergy(obj)
            % Compute turbulent kinetic energy (TKE) amplification ratio from the data
            %
            % Args:
            %     obj (HTRDataAverage): HTRDataAverage object
            %
            % Returns:
            %     TKE_f (float): Turbulent kinetic energy (TKE)
            %     R11_f (float): Reynolds stress R11
            %     R22_f (float): Reynolds stress R22
            %     R33_f (float): Reynolds stress R33
            %     RTT_f (float): Reynolds stress RTT (mean of R22 and R33)

            % Definitions
            rho_avg_ndim = obj.averages.rho_avg ./ obj.averages.rho_avg(obj.indexPreShock);

            % Compute turbulent kinetic energy (TKE)
            obj.averages.TKE = sum(obj.averages.velocity_rms) ./ sum(obj.averages.velocity_rms(:, obj.indexPreShock));
            obj.averages.TKE_f = sum(obj.averages.velocity_frms) ./ sum(obj.averages.velocity_frms(:, obj.indexPreShock)) ./ rho_avg_ndim;

            % Compute Reynolds stress R11
            obj.averages.R11 = obj.averages.velocity_rms(1, :) ./ obj.averages.velocity_rms(1, obj.indexPreShock);
            obj.averages.R11_f = obj.averages.velocity_frms(1, :) ./ obj.averages.velocity_frms(1, obj.indexPreShock) ./ rho_avg_ndim;
            obj.averages.R11_tke = obj.averages.velocity_rms(1, :) ./ (1/3 * sum(obj.averages.velocity_rms));
            obj.averages.R11_ftke = obj.averages.velocity_frms(1, :) ./ (1/3 * sum(obj.averages.velocity_frms)) ./ rho_avg_ndim;

            % Compute Reynolds stress R22
            obj.averages.R22 = obj.averages.velocity_rms(2, :) ./ obj.averages.velocity_rms(2, obj.indexPreShock);
            obj.averages.R22_f = obj.averages.velocity_frms(2, :) ./ obj.averages.velocity_frms(2, obj.indexPreShock) ./ rho_avg_ndim;
            obj.averages.R22_tke = obj.averages.velocity_rms(2, :) ./ (1/3 * sum(obj.averages.velocity_rms));
            obj.averages.R22_ftke = obj.averages.velocity_frms(2, :) ./ (1/3 * sum(obj.averages.velocity_frms)) ./ rho_avg_ndim;

            % Compute Reynolds stress R33
            obj.averages.R33 = obj.averages.velocity_rms(3, :) ./ obj.averages.velocity_rms(3, obj.indexPreShock);
            obj.averages.R33_f = obj.averages.velocity_frms(3, :) ./ obj.averages.velocity_frms(3, obj.indexPreShock) ./ rho_avg_ndim;
            obj.averages.R33_tke = obj.averages.velocity_rms(3, :) ./ (1/3 * sum(obj.averages.velocity_rms));
            obj.averages.R33_ftke = obj.averages.velocity_frms(3, :) ./ (1/3 * sum(obj.averages.velocity_frms)) ./ rho_avg_ndim;

            % Compute Reynolds stress RTT (mean of R22 and R33)
            obj.averages.RTT = sum(obj.averages.velocity_rms(2:3, :)) ./ sum(obj.averages.velocity_rms(2:3, obj.indexPreShock));
            obj.averages.RTT_f = sum(obj.averages.velocity_frms(2:3, :)) ./ sum(obj.averages.velocity_frms(2:3, obj.indexPreShock)) ./ rho_avg_ndim;
            obj.averages.RTT_tke = sum(obj.averages.velocity_rms(2:3, :)) ./ (1/3 * sum(obj.averages.velocity_rms));
            obj.averages.RTT_ftke = sum(obj.averages.velocity_frms(2:3, :)) ./ (1/3 * sum(obj.averages.velocity_frms)) ./ rho_avg_ndim;

            % Assign properties 
            TKE_f = obj.averages.TKE_f;
            R11_f = obj.averages.R11_f;
            R22_f = obj.averages.R22_f;
            R33_f = obj.averages.R33_f;
            RTT_f = obj.averages.RTT_f;
        end

        function indexFarfield = getFarfieldLocation(obj)
            % Get far-field location from the data
            %
            % Args:
            %     obj (HTRDataAverage): HTRDataAverage object
            %
            % Returns:
            %     indexFarfield (float): Index of the far-field location

            % Definitions
            [~, indexFarfield] = max(obj.averages.R11(obj.indexPostShock:end));

            % Apply offset
            indexFarfield = indexFarfield + obj.indexPostShock;

            % Assign properties
            obj.indexFarfield = indexFarfield;
        end

        function variances = getVariances(obj)
            % Compute one-point statistics of the data
            %
            % Note: The correlations are normalized by the pre-shock values
            %
            % Args:
            %     obj (HTRDataAverage): HTRDataAverage object
            %
            % Returns:
            %     variances (float): Variances of the data

            % Compute variances
            obj.averages.rho_variance = obj.averages.rho_rms / obj.averages.rho_avg(obj.indexPreShock)^2;
            obj.averages.pressure_variance =  obj.averages.pressure_frms / (obj.averages.rho_avg(obj.indexPreShock) * obj.averages.SoS_avg(obj.indexPreShock)^2)^2;
            obj.averages.temperature_variance =  obj.averages.temperature_frms / obj.averages.temperature_favg(obj.indexPreShock)^2;
            obj.averages.u_variance = obj.averages.velocity_rms(1, :) / obj.averages.SoS_avg(obj.indexPreShock)^2;
            obj.averages.v_variance = obj.averages.velocity_rms(2, :) / obj.averages.SoS_avg(obj.indexPreShock)^2;
            obj.averages.w_variance = obj.averages.velocity_rms(3, :) / obj.averages.SoS_avg(obj.indexPreShock)^2;

            % Assign properties
            variances.rho = obj.averages.rho_variance;
            variances.pressure = obj.averages.pressure_variance;
            variances.temperature = obj.averages.temperature_variance;
            variances.u = obj.averages.u_variance;
            variances.v = obj.averages.v_variance;
            variances.w = obj.averages.w_variance;
        end

        function correlations = getCorrelations(obj)
            % Compute two-point statistics of the data (correlations)
            % 
            % Note: The correlations are normalized by the pre-shock values
            %
            % Args:
            %     obj (HTRDataAverage): HTRDataAverage object
            %
            % Returns:
            %     correlations (float): Correlations of the data

            % Compute correlations (normalized by the pre-shock values)
            obj.averages.corr_rhoT = (obj.averages.temperature_favg .* obj.averages.rho_avg - obj.averages.rho_avg .* obj.averages.temperature_avg) ./ (obj.averages.rho_avg(obj.indexPreShock) .* obj.averages.temperature_avg(obj.indexPreShock));
            obj.averages.corr_rhoP = (obj.averages.pressure_favg .* obj.averages.rho_avg - obj.averages.rho_avg .* obj.averages.pressure_avg) ./ (obj.averages.rho_avg(obj.indexPreShock) .* obj.averages.SoS_avg(obj.indexPreShock))^2;
            obj.averages.corr_rhoU = (obj.averages.velocity_favg(1, :) .* obj.averages.rho_avg - obj.averages.rho_avg .* obj.averages.velocity_avg(1, :))./ (obj.averages.rho_avg(obj.indexPreShock) .* obj.averages.SoS_avg(obj.indexPreShock));
            obj.averages.corr_rhoV = (obj.averages.velocity_favg(2, :) .* obj.averages.rho_avg - obj.averages.rho_avg .* obj.averages.velocity_avg(2, :)) ./ (obj.averages.rho_avg(obj.indexPreShock) .* obj.averages.SoS_avg(obj.indexPreShock));
            obj.averages.corr_rhoW = (obj.averages.velocity_favg(3, :) .* obj.averages.rho_avg - obj.averages.rho_avg .* obj.averages.velocity_avg(3, :)) ./ (obj.averages.rho_avg(obj.indexPreShock) .* obj.averages.SoS_avg(obj.indexPreShock));
            obj.averages.corr_uT = (obj.averages.uT_avg(1, :) - obj.averages.velocity_avg(1, :) .* obj.averages.temperature_avg) ./ (obj.averages.SoS_avg(obj.indexPreShock) .* obj.averages.temperature_avg(obj.indexPreShock));
            obj.averages.corr_vT = (obj.averages.uT_avg(2, :) - obj.averages.velocity_avg(2, :) .* obj.averages.temperature_avg) ./ (obj.averages.SoS_avg(obj.indexPreShock) .* obj.averages.temperature_avg(obj.indexPreShock));
            obj.averages.corr_wT = (obj.averages.uT_avg(3, :) - obj.averages.velocity_avg(3, :) .* obj.averages.temperature_avg) ./ (obj.averages.SoS_avg(obj.indexPreShock) .* obj.averages.temperature_avg(obj.indexPreShock));

            % Assign properties
            correlations.rhoT = obj.averages.corr_rhoT;
            correlations.rhoP = obj.averages.corr_rhoP;
            correlations.rhoU = obj.averages.corr_rhoU;
            correlations.rhoV = obj.averages.corr_rhoV;
            correlations.rhoW = obj.averages.corr_rhoW;
            correlations.uT = obj.averages.corr_uT;
            correlations.vT = obj.averages.corr_vT;
            correlations.wT = obj.averages.corr_wT;
        end

        function print(obj)
            % Print summary of the average data
            %
            % Args:
            %     obj (HTRDataAverage): HTRDataAverage object
            %
            % Example:
            %     print(obj);

            % Print summary
            fprintf('*********************************************\n');
            fprintf('Summary of the average data\n');
            fprintf('---------------------------------------------\n');
            fprintf('    Parameter              |   Value\n');
            fprintf('---------------------------------------------\n');
            fprintf(' Axis                      |     %s\n', obj.axis);
            fprintf(' Pre-shock location        |   %.4f\n', obj.averages.x(obj.indexPreShock));
            fprintf(' Post-shock location       |   %.4f\n', obj.averages.x(obj.indexPostShock));
            fprintf(' Far-field location        |   %.4f\n', obj.averages.x(obj.indexFarfield));
            fprintf('---------------------------------------------\n');
            fprintf(' Mach number               |   %.4f\n', obj.averages.Ma(obj.indexPreShock));
            fprintf(' Turbulent Mach number     |   %.4f\n', obj.averages.Mt(obj.indexPreShock));
            fprintf(' Kolmogorov length         |   %.4f\n', obj.averages.lengthKolmogorov(obj.indexPreShock));
            fprintf(' Min ratio grid/Kolmogorov |   %.4f\n', min(obj.averages.ratioGridLengthKolmogorov));
            fprintf(' Max ratio grid/Kolmogorov |   %.4f\n', max(obj.averages.ratioGridLengthKolmogorov));
            fprintf(' Dissipation rate          |   %.4f\n', obj.averages.dissipation(obj.indexFarfield));
            fprintf('---------------------------------------------\n');
            fprintf(' TKE amplification ratio   |   %.4f\n', obj.averages.TKE_f(obj.indexFarfield));
            fprintf(' R11 amplification ratio   |   %.4f\n', obj.averages.R11_f(obj.indexFarfield));
            fprintf(' R22 amplification ratio   |   %.4f\n', obj.averages.R22_f(obj.indexFarfield));
            fprintf(' R33 amplification ratio   |   %.4f\n', obj.averages.R33_f(obj.indexFarfield));
            fprintf(' RTT amplification ratio   |   %.4f\n', obj.averages.RTT_f(obj.indexFarfield));
            fprintf('*********************************************\n');
        end

        function ax = plot(obj, xField, yField, varargin)
            % Plot the average data
            %
            % Args:
            %     obj (HTRDataAverage): HTRDataAverage object
            %     xField (char): Name of the x-coordinate property
            %     yField (char): Name of the y-coordinate property
            %
            % Additional key-value arguments:
            %     ax (matlab.graphics.axis.Axes): Axes object to plot the data
            %
            % Returns:
            %     ax (matlab.graphics.axis.Axes): Axes object with the plot

            % Import packages
            import combustiontoolbox.utils.display.*

            % Parse input arguments
            p = inputParser;
            addParameter(p, 'ax', [], @(x) isa(x, 'matlab.graphics.axis.Axes') || isa(x, 'matlab.graphics.layout.TiledChartLayout'));
            addParameter(p, 'xComponent', [], @(x) isnumeric(x) && x > 0 && x <= 3);
            addParameter(p, 'yComponent', [], @(x) isnumeric(x) && x > 0 && x <= 3);
            parse(p, varargin{:});

            % Assign properties
            ax = p.Results.ax;
            xComponent = p.Results.xComponent;
            yComponent = p.Results.yComponent;

            % Definitions
            xValue = obj.averages.(xField);
            yValue = obj.averages.(yField);

            % Check if xComponent and yComponent are provided
            if ~isempty(xComponent)
                xValue = xValue(xComponent, :);
            end

            if ~isempty(yComponent)
                yValue = yValue(yComponent, :);
            end

            % Check if xField and yField are the same
            numCasesX = length(xValue);
            numCasesY = length(yValue);

            if numCasesX < numCasesY
                numCasesY = numCasesX;
            else
                numCasesX = numCasesY;
            end

            xValue = xValue(1:numCasesX);
            yValue = yValue(1:numCasesY);
            
            % Labels
            if obj.FLAG_K0 && strcmp(xField, 'x') && obj.FLAG_SHIFT
                xField = 'k_0 (x - x_{\rm shock})';
            elseif obj.FLAG_K0 && strcmp(xField, 'x')
                xField = 'k_0 x';
            end

            if obj.FLAG_K0 && strcmp(yField, 'x') && obj.FLAG_SHIFT
                yField = 'k_0 (x - x_{\rm shock})';
            elseif obj.FLAG_K0 && strcmp(yField, 'x')
                yField = 'k_0 x';
            end

            % Check if ax is a TiledChartLayout object
            if isa(ax, 'matlab.graphics.layout.TiledChartLayout')
                nexttile(ax); ax = gca; setFigure(ax, obj.plotConfig);
            end

            % Check if ax is provided
            if isempty(ax)
                ax = setFigure(obj.plotConfig);
            else
                
                xLabel = strrep(ax.XLabel.String, '$', '');
                xLabel = strrep(xLabel, ' ', '');

                if strcmp(xLabel, strrep(xField, ' ', ''))
                    yField = '\rm Multiple\ properties';
                end
                
                yLabel = strrep(ax.YLabel.String, '$', '');
                yLabel = strrep(yLabel, ' ', '');

                if strcmp(yLabel, strrep(yField, ' ', ''))
                    xField = '\rm Multiple\ properties';
                end

            end

            % Plot data
            plotFigure(xField, xValue, yField, yValue, 'color', 'auto', 'ax', ax);
            plot(ax, xValue(obj.indexPreShock), yValue(obj.indexPreShock), 'kd', 'MarkerSize', 10, 'linewidth', obj.plotConfig.linewidth);
            plot(ax, xValue(obj.indexPostShock), yValue(obj.indexPostShock), 'kd', 'MarkerSize', 10, 'linewidth', obj.plotConfig.linewidth);
            plot(ax, xValue(obj.indexFarfield), yValue(obj.indexFarfield), 'kd', 'MarkerSize', 10, 'linewidth', obj.plotConfig.linewidth);
            
            ylim(ax, [min(ax.YLim(1), obj.plotYaxisScale(1) * min(min(yValue))), obj.plotYaxisScale(2) * max([yValue(obj.indexPreShock), yValue(obj.indexPostShock), yValue(obj.indexFarfield)])])
        end

    end

    methods (Access = private)

        function applyWeight(obj, property)
            % Apply weight to the raw data
            
            % Check if is numeric
            if ~isnumeric(obj.rawData.(property))
                return
            end

            % Definitions
            numCases = length(obj.rawData.(property)(:, 1));

            % Apply weight
            for i = numCases:-1:1
                obj.averages.(property)(i, :) = obj.rawData.(property)(i, :) ./ obj.weight;
            end

        end

        function rms = getRMS(obj, property)
            % Compute Reynolds root-mean-square (rms) of the data

            % Definitions
            numCases = length(obj.rawData.(property)(:, 1));
            property_avg = strrep(property, '_rms', '_avg');

            % Check if property_avg exists
            if ~isfield(obj.averages, property_avg)
                property = strrep(property, '_rms', '');

                for i = numCases:-1:1
                    rms(i, :) = obj.averages.(property)(i, :).^2;
                end
                
                return
            end

            % Compute Reynolds rms
            for i = numCases:-1:1
                rms(i, :) = obj.averages.(property)(i, :) - obj.averages.(property_avg)(i, :).^2;
            end

        end

        function frms = getFavreRMS(obj, property)
            % Compute Favre root-mean-square (rms) of the data

            % Definitions
            numCases = length(obj.rawData.(property)(:, 1));
            property_favg = strrep(property, '_frms', '_favg');
            
            % Compute Favre rms
            for i = numCases:-1:1
                frms(i, :) = obj.averages.(property)(i, :) - obj.averages.(property_favg)(i, :).^2 ./ obj.averages.rho_avg;
            end
            
        end

    end
    
end
