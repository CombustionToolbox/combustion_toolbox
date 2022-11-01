function [x, y] = smooth_data(x, y, start_point)
    % Prepare data for curve fitting
    [xData, yData] = prepareCurveData(x, y);
    % Set up fittype and options
    ft = fittype('fourier2');
    opts = fitoptions('Method', 'NonlinearLeastSquares');
    opts.Display = 'Off';
    opts.StartPoint = start_point;
    % Fit model to data
    fitresult = fit(xData, yData, ft);
    % Get smooth values
    y = fitresult(xData)';
    x = xData';
    % See results
    figure;
    plot(fitresult, xData, yData);
end