function [tMean, tArray] = timeFunction(f, nFrec)
    % Function to compute mean evaluation time of a given function f
    %
    % Args:
    %      * f (function): Function to be evaluated
    %      * nFrec (float): Number of evaluations
    %
    % Returns:
    %      * tMean (float): Average evaluation time
    %      * tArray (float): Array with the evaluation times
    %
    % Note: The function needs at least one input
    %  
    % Example:
    %      Tuple containing
    %
    %      * [tMean, tArray] = timeFunction(@Example_TP_scoggins2015, 10)

    for i = nFrec:-1:1
        t = tic;
        f(1);
        tArray(i) = toc(t);
    end
    
    tMean = mean(tArray);
end
