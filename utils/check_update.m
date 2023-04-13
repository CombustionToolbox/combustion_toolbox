function [FLAG_UPDATE, message] = check_update(varargin)
    % Check if there is a new release of the Combustion Toolbox
    %
    % Optional Args:
    %    fig (object): UIFigure class
    %
    % Returns:
    %    Tuple containing
    %
    %    * FLAG_UPDATE (bool): FLAG indicating true (false) if there is (not) an update of the Combustion Toolbox
    %    * message (char): Message displayed 
    %
    % Examples:
    %    * [FLAG_UPDATE, message] = check_update();
    %    * [FLAG_UPDATE, message] = check_update(UIFigure);

    % Definitions
    user = 'AlbertoCuadra';
    repo_name = 'combustion_toolbox';

    % Unpack additional inputs
    fig = unpack(varargin{:});

    % Initialization
    FLAG_UPDATE = false;

    % Get current version of CT
    [tag_current, date_current] = get_combustion_toolbox_version();

    % Get latest version of CT on Github
    [tag_repo, git_data] = get_latest_version_github(user, repo_name);
    date_repo = git_data.published_at(1:end - 1);

    % Convert tag to id
    tag_id_current = get_tag_id(tag_current);
    tag_id_repo = get_tag_id(tag_repo);

    % Compare versions
    l_current = length(tag_id_current);
    l_report = length(tag_id_repo);
    N = min(l_current, l_report);
    i = 0;
    
    % Format date
    date_current = datetime(date_current, 'Format', 'dd MMM yyyy');
    date_repo = datetime(date_repo, 'Format', 'uuuu-MM-dd''T''HH:mmXXX', 'TimeZone', 'UTC');
    date_repo.Format = 'dd MMM yyyy';

    % Check latest version
    while ~FLAG_UPDATE && i < N
        i = i + 1;

        if tag_id_current(i) < tag_id_repo(i)
            FLAG_UPDATE = true;
        end

    end

    % Get display message
    [message, title, icon] = get_message(tag_current, tag_repo, date_repo, FLAG_UPDATE);

    % Display modal dialog
    if isa(fig, 'matlab.ui.Figure')
        uialert(fig, message, title, 'Icon', icon);
    else
        S.Interpreter = 'tex';
        S.WindowStyle = 'modal';

        if strcmpi(icon, 'warning')
            icon = 'warn';
            msgbox(['\fontsize{10} ', message], title, icon, S);
            return
        end
        
        icon = 'custom';
        [icondata, iconcmap] = imread('icon_success.ico');
        msgbox(['\fontsize{10} ', message], title, icon, icondata, iconcmap, S); 
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

function [message, title, icon] = get_message(tag_1, tag_2, date_2, FLAG_UPDATE)
    % Get display message
    if FLAG_UPDATE
        msg_1 = 'There is a new release of Combustion Toolbox.';
        title = 'Warning';
        icon = 'warning';
    else
        msg_1 = 'Combustion Toolbox is up to date.';
        title = 'Success';
        icon = 'success';
    end

    message = sprintf(['%s\n\n',...
                       'Current release: %s\n',...
                       'Latest relase: %s\n',...
                       'Date: %s'], msg_1, tag_1, tag_2, date_2);

    % Print message in the command window
    disp(message);
end