function debug_plot_error(it, STOP, lambda)
    % Debug function that plots the error per iteration along with the
    % value of the correction factor

    figure(1); hold on; set(gca,'yscale','log')
    yyaxis left
    plot(1:1:it, STOP, '-', 'Color', '#0072BD');
    yyaxis right
    plot(1:1:it, lambda, '-', 'Color', '#D95319');
end