function position = getMonitorPositionsMATLAB(varargin)
    % Routine that gets the position in pixels of the monitor(s) connected
    % to the device using MATLAB's routines. If no monitor is specified,
    % the position of the main monitor is returned.
    %
    % Optional Args:
    %    monitor_id (float): Get position in pixels for the given monitor
    %
    % Returns:
    %    position (float): Position of the monitor(s)
    %
    % Examples:
    %     * position = get_monitor_positions() returns the position of the main monitor
    %     * position = getMonitorPositions(2) returns the position of the second monitor
    %     * position = getMonitorPositions(1) returns the position of the first monitor
    %     * position = getMonitorPositions(monitor_id) returns the position in pixels of the given monitor

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