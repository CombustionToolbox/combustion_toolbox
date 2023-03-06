function debug_plot_error(it, STOP, varargin)
    % Debug function that plots the error per iteration along with the
    % value of the correction factor

    figure(1); hold on; set(gca, 'yscale', 'log')
    xlabel('Iterations', 'Interpreter', 'latex')
    yyaxis left
    plot(1:1:it, STOP, '-', 'Color', '#0072BD');
    ylabel('Error', 'Interpreter', 'latex')
    
    if nargin > 2
        Delta = varargin{1};
        yyaxis right
        plot(1:1:it, Delta, '-', 'Color', '#D95319');
        ylabel('$\Delta$', 'Interpreter', 'latex')
    end

end
