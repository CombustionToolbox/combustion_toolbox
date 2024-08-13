function [tMean, tArray] = timeFunction(f, nFrec, varargin)
    % Function to compute mean evaluation time of a given function f
    %
    % Args:
    %      f (function): Function to be evaluated
    %      nFrec (float): Number of evaluations
    %
    % Returns:
    %      Tuple containing
    %
    %      * tMean (float): Average evaluation time
    %      * tArray (float): Array with the evaluation times
    %
    % Note: The function needs at least one input
    %  
    % Example:
    %      [tMean, tArray] = timeFunction(@Example_TP_scoggins2015, 10)
    
    % Import packages
    import combustiontoolbox.utils.clearCache

    for i = nFrec:-1:1
        clearCache();
        t = tic;
        f(varargin{:});
        tArray(i) = toc(t);
    end
    
    tMean = mean(tArray);
end
