function FLAG_UPDATE = check_update(varargin)
    % Check if there is a new release of Combustion Toolbox
    %
    % Optional Args:
    %    fig (uifigure): UIFigure class
    %
    % Returns:
    %    FLAG_UPDATE (bool): FLAG indicating true (false) if there is (not) an update of Combustion Toolbox

    % Unpack additional inputs
    fig = unpack(varargin{:});
    % Initialization
    FLAG_UPDATE = true;
    % Get current version of CT
    [tag_1, date_1] = get_combustion_toolbox_version();
    % Get latest version of CT on Github
    [tag_2, git_data] = get_latest_version_github('AlbertoCuadra', 'combustion_toolbox');
    date_2 = git_data.published_at(1:end - 1);
    % Convert tag to id
    tag_id_1 = get_tag_id(tag_1);
    tag_id_2 = get_tag_id(tag_2);
    % Compare versions
    l_1 = length(tag_id_1);
    l_2 = length(tag_id_2);
    N = min(l_1, l_2);
    i = 0;

    while FLAG_UPDATE && i < N
        i = i + 1;

        if tag_id_1(i) > tag_id_2(i)
            FLAG_UPDATE = false;
        end

    end

    % Get display message
    if FLAG_UPDATE
        % Combustion Toolbox is NOT up-to-date
        date_2 = datetime(date_2, 'Format', 'uuuu-MM-dd''T''HH:mmXXX', 'TimeZone', 'UTC');
        date_2.Format = 'dd MMM yyyy';
        message = sprintf(['There is a new release of Combustion Toolbox.\n\n', ...
                            'Latest relase: %s\n', ...
                        'Date: %s'], tag_2, date_2);
        title = 'Warning';
        icon = 'warning';
    else
        % Combustion Toolbox is up-to-date
        date_1 = datetime(date_1, 'Format', 'dd MMM yyyy');
        message = sprintf(['Combustion Toolbox is up to date.\n\n', ...
                            'Latest relase: %s\n', ...
                        'Date: %s'], tag_1, date_1);
        title = 'Success';
        icon = 'success';
    end

    % Display modal dialog
    if isa(fig, 'matlab.ui.Figure')
        uialert(fig, message, title, 'Icon', icon);
    else

        if strcmpi(icon, warning)
            icon = 'warn';
        else
            icon = 'none';
        end

        msgbox(message, title, icon);
    end

end

% SUB-PASS FUNCTIONS
function fig = unpack(varargin)

    if nargin
        fig = varargin{1};
    else
        fig = [];
    end

end

function value = get_tag_id(tag)
    ind = regexp(tag, '\d');
    value = tag(ind);
end
