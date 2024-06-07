function position = getMonitorPositions(varargin)
    % Routine that gets the position in pixels of the monitor(s) connected
    % to the device using Java (default) or MATLAB's routines. If no monitor
    % is specified, the position of the main monitor is returned.
    %
    % Optional Args:
    %     monitor_id (float): Get position in pixels for the given monitor
    %
    % Returns:
    %     position (float): Position of the monitor(s)
    %
    % Examples:
    %     * position = get_monitor_positions() returns the position of the main monitor
    %     * position = getMonitorPositions(2) returns the position of the second monitor
    %     * position = getMonitorPositions(1) returns the position of the first monitor
    %     * position = getMonitorPositions(monitor_id) returns the position in pixels of the given monitor
    %
    % Notes:
    %     This function first tries using Java to get the screen size values that do 
    %     not take into account ui scaling, so it will return the proper screen size 
    %     values. Otherwise, it gets the screen position using MATLAB's routines
    
    % Get the position of the monitor using Java
    try
        dimensions = java.awt.Toolkit.getDefaultToolkit().getScreenSize();
        position_monitors = [1, 1, dimensions.width, dimensions.height];
        if nargin
            % Get position for the given monitor
            try
                monitor_id = varargin{1};
                % Get local graphics environment
                env = java.awt.GraphicsEnvironment.getLocalGraphicsEnvironment();
                % Get screen device objects
                devices = env.getScreenDevices();
                % Get initial points
                position = get_monitor_positions_MATLAB(varargin{:});
                % Get screen size for the selected monitor
                position = [position(1), position(2), ...
                    devices(monitor_id).getDisplayMode().getWidth(),...
                    devices(monitor_id).getDisplayMode().getHeight()];
            catch
                monitor_id = 1;
                position = position_monitors(monitor_id, :);
            end

        else
            % Find main monitor
            position = position_monitors;
        end
    
    % Get the position of the monitor using MATLAB's routines
    catch
        position = getMonitorPositionsMATLAB(varargin{:});
    end

end