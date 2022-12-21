function [x, y] = smooth_data(x, y, varargin)
    % Smooth data using Fourier NonlinearLeastSquares method
    % 
    % Args:
    %     x (float): data in the x direction
    %     y (float): data in the y direction
    %
    % Optional Args:
    %     start_point (float): initial point of the fit
    %
    % Return:
    %     Tuple
    %     * x (float): smooth data in the x direction
    %     * y (float): smooth data in the y direction
    
    % Definitions
    start_point = [0 0 0 0 0 0];
    method = 'fourier2';

    % Unpack
    if nargin > 2
        start_point = varargin{1};
    end
    
    if nargin > 3
        method = varargin{2};
    end
    
    % Prepare data for curve fitting
    [xData, yData] = prepareCurveData(x, y);
    % Set up fittype and options
    ft = fittype(method);
    opts = fitoptions('Method', 'NonlinearLeastSquares');
    opts.Display = 'Off';
    opts.StartPoint = start_point;
    % Fit model to data
    fitresult = fit(xData, yData, ft, opts);
    % Get smooth values
    y = fitresult(xData)';
    x = xData';
    % See results
    figure;
    plot(fitresult, xData, yData);
end