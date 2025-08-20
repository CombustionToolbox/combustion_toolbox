function x = clusteredMesh1D(xDense, xCoarse, nDense, nCoarse, varargin)
    % Create a 1D nonuniform mesh with dense and coarse regions.
    %
    % Args:
    %    xDense (float): Dense region [xDense(1), xDense(2)]
    %    xCoarse (float): Coarse region [xCoarse(1), xCoarse(2)]
    %    nDense (float): Number of points in the dense region
    %    nCoarse (float): Number of points in the coarse region
    %
    % Optional key-value pairs:
    %    * 'SpacingDense'  : 'log' (default) or 'linear'
    %    * 'SpacingCoarse' : 'log' (default) or 'linear'
    %    * 'Sort'          : true (default), sorts final mesh
    % 
    % Returns:
    %    x (float): 1D array of nonuniform mesh points
    %
    % Examples:
    %    * x = clusteredMesh1D([1, 1.2], [1.2, 10], 50, 50);
    %    * x = clusteredMesh1D([1, 1.2], [1.2, 10], 50, 50, 'SpacingDense', 'log', 'SpacingCoarse', 'linear');

    % Parse optional args
    p = inputParser;
    addParameter(p, 'SpacingDense', 'log');
    addParameter(p, 'SpacingCoarse', 'log');
    addParameter(p, 'Sort', true);
    parse(p, varargin{:});
    opt = p.Results;

    % Get dense region
    if strcmpi(opt.SpacingDense, 'log')

        if any(xDense <= 0)
            error('Log spacing requires xDense > 0.');
        end

        x1 = logspace(log10(xDense(1)), log10(xDense(2)), nDense);
    else
        x1 = linspace(xDense(1), xDense(2), nDense);
    end

    % Get coarse region
    if strcmpi(opt.SpacingCoarse, 'log')

        if any(xCoarse <= 0)
            error('Log spacing requires xCoarse > 0.');
        end

        x2 = logspace(log10(xCoarse(1)), log10(xCoarse(2)), nCoarse);
    else
        x2 = linspace(xCoarse(1), xCoarse(2), nCoarse);
    end

    % Get unique points
    x = unique([x1, x2]);

    % Optionally sort
    if opt.Sort
        x = sort(x);
    end

end
