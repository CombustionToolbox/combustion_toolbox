classdef TurbulenceSpectra < handle
    % The :mat:func:`TurbulenceSpectra` class provides methods for computing
    % turbulence spectra. Supports spherical and cross-plane averaging.

    properties
        averaging = 'spherical'  % Type of averaging ('spherical', 'spherical2D', 'crossplane', 'crossplane2D')
        axis = 'x'               % Axis for cross-plane averaging ('x', 'y', 'z')
        fps = 60                 % Frame rate for movie (temporal will be added into PlotConfig)
        time                     % Elapsed time
        plotConfig               % PlotConfig object
        FLAG_TIME = true         % Flag to print elapsed time
        FLAG_REPORT = false      % Flag to print predefined plots
    end

    methods

        function obj = TurbulenceSpectra(varargin)
            % Constructor for TurbulenceSpectra
            %
            % Optional Args:
            %     averaging (char): Type of averaging ('spherical', 'spherical2D', 'crossplane', 'crossplane2D')
            %     axis (char): Axis for cross-plane averaging ('x', 'y', 'z')
            %
            % Example:
            %     analyzer = TurbulenceSpectra('averaging', 'crossplane', 'axis', 'z');

            % Default
            defaultPlotConfig = combustiontoolbox.utils.display.PlotConfig();
            defaultPlotConfig.xscale = 'log'; defaultPlotConfig.yscale = 'log';
            defaultPlotConfig.innerposition = [0.15 0.15 0.35 0.5];
            defaultPlotConfig.outerposition = [0.15 0.15 0.35 0.5];

            % Parse input arguments
            p = inputParser;
            addParameter(p, 'averaging', obj.averaging, @(x) ischar(x) && ismember(lower(x), {'spherical', 'spherical2d', 'crossplane', 'crossplane2d'}));
            addParameter(p, 'axis', obj.axis, @(x) ischar(x) && ismember(lower(x), {'x', 'y', 'z'}));
            addParameter(p, 'plotConfig', defaultPlotConfig, @(x) isa(x, 'combustiontoolbox.utils.display.PlotConfig'));
            addParameter(p, 'FLAG_TIME', obj.FLAG_TIME, @islogical);
            addParameter(p, 'FLAG_REPORT', obj.FLAG_REPORT, @islogical);
            parse(p, varargin{:});

            obj.averaging = lower(p.Results.averaging);
            obj.axis = lower(p.Results.axis);
            obj.plotConfig = p.Results.plotConfig;
            obj.FLAG_TIME = p.Results.FLAG_TIME;
            obj.FLAG_REPORT = p.Results.FLAG_REPORT;
        end

        function obj = set(obj, property, value, varargin)
            % Set properties of the TurbulenceSpectra object
            %
            % Args:
            %     obj (TurbulenceSpectra): TurbulenceSpectra object
            %     property (char): Property name
            %     value (float): Property value
            %
            % Optional Args:
            %     * property (char): Property name
            %     * value (float): Property value
            %
            % Returns:
            %     obj (TurbulenceSpectra): TurbulenceSpectra object with updated properties
            %
            % Examples:
            %     * set(TurbulenceSpectra(), 'averaging', 'crossplane');
            %     * set(TurbulenceSpectra(), 'averaging', 'crossplane', 'axis', 'z');
            %     * set(TurbulenceSpectra(), 'averaging', 'crossplane', 'axis', 'z', 'FLAG_TIME', false);
            
            varargin = [{property, value}, varargin{:}];

            for i = 1:2:length(varargin)
                % Assert that the property exists
                assert(isprop(obj, varargin{i}), 'Property not found');

                % Set property
                obj.(varargin{i}) = varargin{i + 1};
            end

        end

        function [EK_avg, k] = getEnergySpectra(obj, fluctuation)
            % Compute the energy spectra of a 3D fluctuation field
            %
            % Args:
            %     obj (TurbulenceSpectra): TurbulenceSpectra object
            %     fluctuation (float): 3D fluctuation field
            %
            % Returns:
            %     EK_avg (float): Averaged energy spectra
            %     k (float): Wavenumber along the homogeneous direction
            %
            % Example:
            %     [EK_avg, k] = getEnergySpectra(TurbulenceSpectra(), fluctuation);

            % Check type of fluctuation
            if isa(fluctuation, 'combustiontoolbox.turbulence.VelocityField')
                [EK_avg, k] = getEnergySpectraVelocity(obj, fluctuation);
                return
            end

            % Check input
            assert(isnumeric(fluctuation), 'Input must be a 3D array.');

            % Start timer
            obj.time = tic;

            % Set the appropriate averaging function
            switch lower(obj.averaging)
                case 'spherical'
                    % Use 3D FFT for spherical averaging
                    F = fftn(fluctuation) / numel(fluctuation);
                    EK = fftshift(abs(F).^2);
                    % Set the function for spherical averaging
                    averagingFunction = @obj.getSphericallyAveragedSpectra;
                case 'spherical2d'
                    % Use 3D FFT for spherical averaging
                    F = fftn(fluctuation) / numel(fluctuation);
                    EK = fftshift(abs(F).^2);
                    % Set the function for spherical averaging
                    averagingFunction = @obj.getSphericallyAveragedSpectra2D;
                case 'crossplane'
                    % Use 2D FFT for cross-plane averaging
                    EK = obj.getCrossplaneFFT(fluctuation, obj.axis);
                    % Set the function for cross-plane averaging
                    averagingFunction = @obj.getCrossplaneAveragedSpectra;
                case 'crossplane2d'
                    % Use 2D FFT for cross-plane averaging
                    EK = obj.getCrossplaneFFT(fluctuation, obj.axis);
                    % Set the function for cross-plane averaging
                    averagingFunction = @obj.getCrossplaneAveragedSpectra2D;
            end
            
            % Compute the energy spectra based on the averaging type
            [EK_avg, k] = averagingFunction(EK);

            % Time elapsed
            obj.time = toc(obj.time);
            
            % Print elapsed time
            printTime(obj);
        end

        function [EK_avg, k] = getEnergySpectraVelocity(obj, velocity)
            % Compute the energy spectrum of a 3D velocity field
            %
            % Args:
            %     obj (TurbulenceSpectra): TurbulenceSpectra object
            %     velocity (VelocityField): Velocity field as a VelocityField object, struct, or 4D matrix
            %
            % Returns:
            %     EK_avg (float): Averaged energy spectra
            %     k (float): Wavenumber along the homogeneous direction
            %
            % Example:
            %     [EK_avg, k] = getEnergySpectraVelocity(TurbulenceSpectra(), velocity);
            
            % Import packages
            import combustiontoolbox.common.Units.convertData2VelocityField

            % Check input
            assert(isa(velocity, 'combustiontoolbox.turbulence.VelocityField'), 'Input must be a VelocityField object.');

            % Start timer
            obj.time = tic;

            % Convert input to VelocityField object
            velocity = convertData2VelocityField(velocity);

            % Set the appropriate averaging function
            switch lower(obj.averaging)
                case 'spherical'
                    % Use 3D FFT for spherical averaging
                    U = fftn(velocity.u) / numel(velocity.u);
                    V = fftn(velocity.v) / numel(velocity.v);
                    W = fftn(velocity.w) / numel(velocity.w);
                    EK = 0.5 * fftshift(abs(U).^2 + abs(V).^2 + abs(W).^2);
                    % Set the function for spherical averaging
                    averagingFunction = @obj.getSphericallyAveragedSpectra;
                case 'spherical2d'
                    % Use 3D FFT for spherical averaging
                    U = fftn(velocity.u) / numel(velocity.u);
                    V = fftn(velocity.v) / numel(velocity.v);
                    W = fftn(velocity.w) / numel(velocity.w);
                    EK = 0.5 * fftshift(abs(U).^2 + abs(V).^2 + abs(W).^2);
                    % Set the function for spherical averaging
                    averagingFunction = @obj.getSphericallyAveragedSpectra2D;
                case 'crossplane'
                    % Use 2D FFT for cross-plane averaging
                    EK = obj.getCrossplaneFFTVelocity(velocity, obj.axis);
                    % Set the function for cross-plane averaging
                    averagingFunction = @obj.getCrossplaneAveragedSpectra;
                case 'crossplane2d'
                    % Use 2D FFT for cross-plane averaging
                    EK = obj.getCrossplaneFFTVelocity(velocity, obj.axis);
                    % Set the function for cross-plane averaging
                    averagingFunction = @obj.getCrossplaneAveragedSpectra2D;
            end
            
            % Compute the energy spectra based on the averaging type
            [EK_avg, k] = averagingFunction(EK);

            % Time elapsed
            obj.time = toc(obj.time);
            
            % Print elapsed time
            printTime(obj);
        end

        function ax = plot(obj, k, EK, varargin)
            % Plot results
            %
            % Args:
            %     obj (TurbluenceSpectra): TurbulenceSpectra object
            %     k (float): Wavenumber along the homogeneous direction
            %     EK (float): Energy spectra
            %
            % Optional Args:
            %     * EK_i (float): Array of energy spectra
            %
            % Returns:
            %     ax (Axes): Axes object
            %
            % Examples:
            %     * ax = plot(TurbulenceSpectra(), k, EK);
            %     * ax = plot(TurbulenceSpectra(), k, EK1, EK2);

            % Import packages
            import combustiontoolbox.utils.display.*
            
            % Definitions
            numCases = length(varargin) + 1;
            numWaves = length(k);

            % Plot spectra
            ax = setFigure(obj.plotConfig);
            plotFigure('k', k, 'E(k)', EK(1:numWaves), 'color', 'auto', 'ax', ax);

            for i = 1:numCases-1
                EK = varargin{i};
                plotFigure('k', k, 'E(k)', EK(1:numWaves), 'color', 'auto', 'ax', ax);
            end

            % legend(ax, {'solenoidal', 'dilatational', 'total'}, 'interpreter', 'latex');
            % plotFigure('k', k, 'E(k)', 0.1 * k.^(-5/3), 'linestyle', '--', 'color', [0 0 0], 'ax', ax);
        end

        function ax = plotPDF(obj, pdf, edges)
            % Plot the probability density function (PDF)
            %
            % Args:
            %     obj (TurbulenceSpectra): TurbulenceSpectra object
            %     pdf (float): Probability density function
            %     edges (float): Bin edges
            %
            % Returns:
            %     ax (Axes): Axes object
            %
            % Example:
            %     ax = plotPDF(obj, pdf, edges);

            % Import packages
            import combustiontoolbox.utils.display.*

            % Definitions
            centers = (edges(1:end-1) + edges(2:end)) / 2;

            % Backup current plot configuration
            originalXScale = obj.plotConfig.xscale;
            originalYScale = obj.plotConfig.yscale;

            % Set plot configuration to linear scales
            obj.plotConfig.xscale = 'linear';
            obj.plotConfig.yscale = 'linear';

            % Plot the PDF
            ax = setFigure(obj.plotConfig);
            plotFigure('x', centers, '\rm{p.d.f.}', pdf, 'color', 'auto', 'ax', ax);

            % Restore original plot configuration
            obj.plotConfig.xscale = originalXScale;
            obj.plotConfig.yscale = originalYScale;
        end

        function ax = plotJointPDF(obj, pdf, edgesX, edgesY)
            % Plot the joint probability density function (PDF)
            %
            % Args:
            %     obj (TurbulenceSpectra): TurbulenceSpectra object
            %     pdf (float): Joint probability density function
            %     edgesX (float): Bin edges for fluctuationX
            %     edgesY (float): Bin edges for fluctuationY
            %
            % Returns:
            %     ax (Axes): Axes object
            %
            % Example:
            %     ax = plotJointPDF(obj, pdf, edgesX, edgesY);

            % Import packages
            import combustiontoolbox.utils.display.*

            % Definitions
            xCenters = (edgesX(1:end-1) + edgesX(2:end)) / 2;
            yCenters = (edgesY(1:end-1) + edgesY(2:end)) / 2;
            [X, Y] = meshgrid(xCenters, yCenters);

            % Backup current plot configuration
            originalXScale = obj.plotConfig.xscale;
            originalYScale = obj.plotConfig.yscale;

            % Set plot configuration to linear scales
            obj.plotConfig.xscale = 'linear';
            obj.plotConfig.yscale = 'linear';

            % Plot the joint PDF
            ax = setFigure(obj.plotConfig);
            contourf(ax, X, Y, log10(pdf)');
            c = colorbar(ax, 'linewidth', obj.plotConfig.linewidth, 'fontsize', obj.plotConfig.fontsize - 4);
            c.Label.Interpreter = 'latex'; c.TickLabelInterpreter = 'latex';
            c.Label.String = '$\log_{10}(\rm{p.d.f.})$';

            % Restore original plot configuration
            obj.plotConfig.xscale = originalXScale;
            obj.plotConfig.yscale = originalYScale;
        end

        function movie(obj, k, EK, varargin)
            % Create a movie by plotting slices of the spectra
            %
            % Args:
            %     obj (TurbulenceSpectra): TurbulenceSpectra object
            %     k (float): Wavenumber along the homogeneous direction
            %     EK (float): Energy spectra (2D array)
            %
            % Example:
            %     obj.movie(k, EK);

            % Import packages
            import combustiontoolbox.utils.display.*

            if ~isnumeric(EK) || size(EK, 1) <= 1
                error('EK must be a 2D array with multiple slices.');
            end

            % Definitions
            numCases = length(varargin) + 1;
            numSlices = size(EK, 1);
            numWaves = length(k);

            % Plot
            ax = setFigure(obj.plotConfig);
            
            % Temporal
            ylim([1e-15, 1e0]);

            % Movie
            for sliceIndex = 1:numSlices
                % Clear previous plot
                cla(ax);

                % Plot current slice
                plotFigure('k', k, 'E(k)', EK(sliceIndex, 1:numWaves), 'color', 'auto', 'ax', ax);

                for i = 1:numCases-1
                    EK_aux = varargin{i};
                    plotFigure('k', k, 'E(k)', EK_aux(sliceIndex, 1:numWaves), 'color', 'auto', 'ax', ax);
                end

                title(ax, sprintf('Slice %d of %d', sliceIndex, numSlices), 'Interpreter', 'latex');

                % Pause for the desired frame rate
                pause(1 / obj.fps);
            end
            
        end

        function printTime(obj)
            % Print execution time
            %
            % Args:
            %     obj (TurbulenceSpectra): TurbulenceSpectra object
            
            if ~obj.FLAG_TIME
                return
            end

            % Definitions
            operationName = 'energy spectra';

            % Print elapsed time
            fprintf('\nElapsed time for %s: %.5f seconds\n', operationName, obj.time);
        end

    end

    methods (Static, Access = public)
        
        function [pdf, edges] = getPDF(fluctuation)
            % Compute the probability density function (PDF) of a 3D fluctuation field
            %
            % Args:
            %     fluctuation (float): 3D fluctuation field
            %
            % Returns:
            %     pdf (float): Probability density function
            %     edges (float): Bin edges
            %
            % Example:
            %     pdf = getPDF(TurbulenceSpectra(), fluctuation);

            % Check input
            assert(isnumeric(fluctuation), 'Input must be a 3D array.');
            
            % Compute the PDF
            [pdf, edges] = histcounts(fluctuation(:), 'Normalization', 'pdf');
        end

        function [pdf, edgesX, edgesY] = getJointPDF(fluctuationX, fluctuationY)
            % Compute the joint probability density function (PDF) of two 3D fluctuation fields
            %
            % Args:
            %     fluctuationX (float): 3D fluctuation field
            %     fluctuationY (float): 3D fluctuation field 
            %
            % Returns:
            %     pdf (float): Joint probability density function
            %     edgesX (float): Bin edges for fluctuationX
            %     edgesY (float): Bin edges for fluctuationY
            %
            % Example:
            %     pdf = getJointPDF(TurbulenceSpectra(), fluctuationX, fluctuationY);

            % Check input
            assert(isnumeric(fluctuationX) && isnumeric(fluctuationY), 'Input must be 3D arrays.');

            % Compute the joint PDF
            [pdf, edgesX, edgesY] = histcounts2(fluctuationX(:), fluctuationY(:), 'Normalization', 'pdf');
        end

        function L = getIntegralLengthScale(EK, k)
            % Compute the integral length scale from the energy spectra
            %
            % Args:
            %     EK (float): Energy spectra
            %     k (float): Wavenumber along the homogeneous direction
            %
            % Returns:
            %     L (float): Integral length scale
            %
            % Example:
            %     L = getIntegralLengthScale(EK, k);

            % Check input
            assert(isnumeric(EK) && isnumeric(k), 'Input must be numeric arrays.');

            % Flatten arrays
            k = k(:);
            EK = EK(:);

            % Remove zero values
            FLAG_REMOVE = k == 0;
            k(FLAG_REMOVE) = [];
            EK(FLAG_REMOVE) = [];

            % Include only wavenumbers contained in the energy spectra
            EK = EK(1:length(k));

            % Compute the integral length scale
            L = 3/4 * pi * trapz(EK ./ k) / trapz(EK);
        end

    end

    methods (Static, Access = private)

        function EK = getCrossplaneFFT(fluctuation, axisType)
            % Compute 2D FFT on homogeneous planes based on axisType
            %
            % Args:
            %   fluctuation (float): 3D fluctuation field
            %   axisType (char): Axis for cross-plane averaging ('x', 'y', 'z')
            %
            % Returns:
            %   EK (float): 3D array representing the 2D energy spectra computed on the
            %               homogeneous planes, sliced along the inhomogeneous axis.
            %               The first two dimensions correspond to the homogeneous plane,
            %               while the third dimension represents the slices along the
            %               inhomogeneous axis.
            %
            % Examples:
            %   * EK = getCrossplaneFFT(fluctuation, 'x');
            %   * EK = getCrossplaneFFT(fluctuation, 'y');
            %   * EK = getCrossplaneFFT(fluctuation, 'z');
            
            % Map axis to slicing dimension
            sliceDim = find('xyz' == axisType);
            dims = 1:3;
            homogeneousDims = setdiff(dims, sliceDim);
        
            % Permute dimensions to bring slice dimension last
            permutedF = permute(fluctuation, [homogeneousDims, sliceDim]);
        
            % Extract slice sizes and initialize EK
            [size1, size2, numSlices] = size(permutedF);
            EK = zeros(size1, size2, numSlices);
        
            % Compute 2D FFT slice by slice
            for i = 1:numSlices
                % Extract slices
                F = fft2(squeeze(permutedF(:, :, i))) / numel(permutedF(:, :, i));
        
                % Compute energy spectra for the slice
                EK(:, :, i) = fftshift(abs(F).^2);
            end
        
        end

        function EK = getCrossplaneFFTVelocity(velocity, axisType)
            % Compute 2D FFT on homogeneous planes based on axisType
            %
            % Args:
            %   u (float): 3D velocity field in the x-direction
            %   v (float): 3D velocity field in the y-direction
            %   w (float): 3D velocity field in the z-direction
            %   axisType (char): Axis for cross-plane averaging ('x', 'y', 'z')
            %
            % Returns:
            %   EK (float): 3D array representing the 2D energy spectra computed on the
            %               homogeneous planes, sliced along the inhomogeneous axis.
            %               The first two dimensions correspond to the homogeneous plane,
            %               while the third dimension represents the slices along the
            %               inhomogeneous axis.
            %
            % Examples:
            %   * EK = getCrossplaneFFTVelocity(velocity, 'x');
            %   * EK = getCrossplaneFFTVelocity(velocity, 'y');
            %   * EK = getCrossplaneFFTVelocity(velocity, 'z');
            
            % Map axis to slicing dimension
            sliceDim = find('xyz' == axisType);
            dims = 1:3;
            homogeneousDims = setdiff(dims, sliceDim);
        
            % Permute dimensions to bring slice dimension last
            permutedU = permute(velocity.u, [homogeneousDims, sliceDim]);
            permutedV = permute(velocity.v, [homogeneousDims, sliceDim]);
            permutedW = permute(velocity.w, [homogeneousDims, sliceDim]);
        
            % Extract slice sizes and initialize EK
            [size1, size2, numSlices] = size(permutedU);
            EK = zeros(size1, size2, numSlices);
        
            % Compute 2D FFT slice by slice
            for i = 1:numSlices
                % Extract slices
                U = fft2(squeeze(permutedU(:, :, i))) / numel(permutedU(:, :, i));
                V = fft2(squeeze(permutedV(:, :, i))) / numel(permutedV(:, :, i));
                W = fft2(squeeze(permutedW(:, :, i))) / numel(permutedW(:, :, i));
        
                % Compute energy spectra for the slice
                EK(:, :, i) = 0.5 * fftshift(abs(U).^2 + abs(V).^2 + abs(W).^2);
            end
        
        end
        
        function [EK_avg, k] = getCrossplaneAveragedSpectra(EK)
            % Compute the cross-plane averaged energy spectra of a 3D field.
            % The homogeneous plane are the first two dimensions of the 3D field.
            %
            % Args:
            %   EK (float): 3D energy spectra (combined from all velocity components)
            %
            % Returns:
            %   EK_avg (float): Radially averaged energy spectra for the slices
            %   k (float): Wavenumber along the homogeneous direction
            %
            % Example:
            %   [EK_avg, k] = getCrossplaneAveragedSpectra(EK);
                
            % Get the center of the homogeneous plane
            [NX, NY, numSlices] = size(EK);
            centerX = floor(NX / 2) + 1;
            centerY = floor(NY / 2) + 1;
            [X, Y] = ndgrid(1:NX, 1:NY);
            
            % Compute radii for radial averaging
            radii = round(sqrt((X  - centerX).^2 + (Y - centerY).^2)) + 1;
            maxRadius = max(radii(:));
            
            % Radially average slices
            EK_avg = zeros(numSlices, maxRadius);
            for i = 1:numSlices
                % Extract slice
                sliceData = EK(:, :, i);
        
                % Radially average
                EK_avg(i, :) = accumarray(radii(:), sliceData(:), [maxRadius, 1], @sum, 1e-50);
            end
        
            % Compute wavenumber
            fftResult = fft(EK(:, 1, 1));
            numFreqs = ceil(length(fftResult) / 2);
            k = sqrt(2) * (0:numFreqs - 1);
        end

        function [EK_avg2D, k] = getCrossplaneAveragedSpectra2D(EK)
            % Compute the cross-plane averaged 2D energy spectra of a 3D field.
            % The homogeneous plane corresponds to the first two dimensions.
            %
            % Args:
            %   EK (float): 3D energy spectra (combined from all velocity components)
            %
            % Returns:
            %   EK_avg2D (float): 2D radially averaged energy spectra
            %   k (float): Wavenumber along homogeneous directions
            %
            % Example:
            %   [EK_avg2D, k] = getCrossplaneAveragedSpectra2D(EK);
        
            % Get the size of the homogeneous plane
            [NX, NY, numSlices] = size(EK);
        
            % Define the center of the homogeneous plane
            centerX = floor(NX / 2) + 1;
            centerY = floor(NY / 2) + 1;
        
            % Create the grid for radial averaging
            [X, Y] = ndgrid(1:NX, 1:NY);
            radiiX = round(abs(X - centerX)) + 1;
            radiiY = round(abs(Y - centerY)) + 1;
        
            % Accumulate energy spectra into 2D grid
            waveIndices = [radiiX(:), radiiY(:)];
            maxWaveIndicesX = max(waveIndices(:, 1));
            maxWaveIndicesY = max(waveIndices(:, 2));
            
            % Radially average slices
            EK_avg2D = zeros(maxWaveIndicesX, maxWaveIndicesY, numSlices);
            for i = 1:numSlices
                % Extract slice
                sliceData = squeeze(EK(:, :, i));

                % Radially average
                EK_avg2D(:, :, i) = accumarray(waveIndices, sliceData(:), ...
                    [maxWaveIndicesX, maxWaveIndicesY], @sum, 1e-50);
            end
        
            % Compute wavenumber
            fftResult = fft(EK(:, 1, 1));
            numFreqs = ceil(length(fftResult) / 2);
            k = sqrt(2) * (0:numFreqs - 1);
        end
        
        function [EK_avg, k] = getSphericallyAveragedSpectra(EK)
            % Compute the spherically averaged energy spectra of a 3D field
            %
            % Args:
            %   EK (float): 3D energy spectra
            %
            % Returns:
            %   EK_avg (float): Spherically averaged energy spectra
            %
            % Example:
            %   EK_avg = getSphericallyAveragedSpectra(EK);
        
            % Get the center of the 3D field
            [NX, NY, NZ] = size(EK);
            centerX = floor(NX / 2) + 1;
            centerY = floor(NY / 2) + 1;
            centerZ = floor(NZ / 2) + 1;
        
            % Compute radii for radial averaging
            [X, Y, Z] = ndgrid(1:NX, 1:NY, 1:NZ);
            radii = round(sqrt((X - centerX).^2 + (Y - centerY).^2 + (Z - centerZ).^2)) + 1;
            maxRadius = max(radii(:));

            % Radially average the energy spectra 
            EK_avg = accumarray(radii(:), EK(:), [maxRadius, 1], @sum, 1e-50);
        
            % Compute wavenumber
            fftResult = fft(EK(:, 1, 1));
            numFreqs = ceil(length(fftResult) / 2);
            k = 0:numFreqs - 1;
        end

        function [EK_avg2D, k] = getSphericallyAveragedSpectra2D(EK)
            % Compute the spherically averaged 2D energy spectra of a 3D field
            %
            % Args:
            %   EK (float): 3D energy spectra
            %
            % Returns:
            %   EK_avg2D (float): Spherically averaged 2D energy spectra
            %   k (float): Wavenumber along the homogeneous direction
            %
            % Example:
            %   [EK_avg2D, k] = getSphericallyAveragedSpectra2D(EK);
            
            % Get the center of the 3D field
            [NX, NY, NZ] = size(EK);
            centerX = floor(NX/2) + 1;
            centerY = floor(NY/2) + 1;

            % Generate grid indices for the entire 3D array
            [X, Y, ~] = ndgrid(1:NX, 1:NY, 1:NZ);

            % Compute wavenumber indices relative to the center
            KX = round(abs(X - centerX));
            KY = round(abs(Y - centerY));

            % Map the 3D energy spectrum to a 2D grid using accumulation
            waveIndices = [KX(:), KY(:)] + 1;

            % Accumulate the energy spectrum in a 2D matrix
            maxWaveIndicesX = max(waveIndices(:, 1));
            maxWaveIndicesY = max(waveIndices(:, 2));
            EK_avg2D = accumarray(waveIndices, EK(:), [maxWaveIndicesX, maxWaveIndicesY], @sum, 1e-50);

            % Compute wavenumber
            fftResult = fft(EK(:, 1, 1));
            numFreqs = ceil(length(fftResult) / 2);
            k = sqrt(2) * (0:numFreqs - 1);
        end

    end

end