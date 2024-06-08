function ax = plot_hugoniot(self, mix1, mix2, varargin)
    % Plot the Hugoniot curve for a a given pre-shock state (mix1)
    % and post-shock state (mix2)
    %
    % Args:
    %     self (struct): Data of the mixture, conditions, and databases
    %     mix1 (struct): Pre-shock mixture
    %     mix2 (struct): Post-shock mixture
    %
    % Optional Args:
    %     ax (object): Axis handle to plot on
    %
    % Returns:
    %     ax (object): Axis handle to plot on
    %
    % Examples:
    %     * ax = plot_hugoniot(self, mix1, mix2)
    %     * ax = plot_hugoniot(self, mix1, mix2, ax)

    
    % Definitions
    self.Misc.config.grid = 'on';
    self.Misc.config.yscale = 'log';
    self.Misc.config.axis_x = [0 , 1];
    self.Misc.config.axis_y = [1, 1e4];
    
    % Unpack
    if nargin > 3
        ax = varargin{1};
        legend
    else
        ax = set_figure(self.Misc.config);
    end

    % Get jump conditions
    rho1 = cell2vector(mix1, 'rho');
    rho2 = cell2vector(mix2, 'rho');
    p1 = cell2vector(mix1, 'p');
    p2 = cell2vector(mix2, 'p');
    R = rho2 ./ rho1;
    P = p2 ./ p1;
    
    % Plot
    ax = plot_figure('$R^{-1}$', 1 ./ R, '$P$', P, 'ax', ax, 'color', 'auto');
end