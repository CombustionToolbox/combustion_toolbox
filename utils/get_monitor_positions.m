function position = get_monitor_positions(varargin)
    % Routine that gets the position in pixels of the monitor(s) connected
    % to the device
    %
    % Optional Args:
    %    monitor_id (float): Get position in pixels for the given monitor
    %
    % Returns:
    %    position (float): Position of the monitor(s)

    
    % Try first using java that does not take into account ui scaling,
    % so it will return the proper screen size values
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
    
    % Otherwise, get screen position using MATLAB's routines
    catch
        position = get_monitor_positions_MATLAB(varargin{:});

    end

end

function position = get_monitor_positions_MATLAB(varargin)
    position_monitors = get(0, 'MonitorPositions');

    if nargin
        % Get position for the given monitor
        try
            monitor_id = varargin{1};
            position = position_monitors(monitor_id, :);
        catch
            monitor_id = 1;
            position = position_monitors(monitor_id, :);
        end

    else
        % Find main monitor
        ind_main = find(position_monitors == 1, 1);
        position = position_monitors(ind_main, :);
    end
end