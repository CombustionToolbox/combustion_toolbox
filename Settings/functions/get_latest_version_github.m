function [release, git_data] = get_latest_version_github(user, repo_name)
    % Get latest version from a repository from Github.
    %
    % Args:
    %    user (char): Username of the owner of the repository
    %    repo_name (char): Name of the repository
    %
    % Returns:
    %    Tuple containing
    %
    %    - release (char): Release tag (latest)
    %    - git_data (struct): Body data of the request

    % Define request
    git_api = 'https://api.github.com/repos/';
    url = [git_api, user, '/', repo_name, '/releases/latest'];
    request = matlab.net.http.RequestMessage;
    uri = matlab.net.URI(url);
    % Send request
    response = send(request, uri);
    % Get data
    git_data = response.Body.Data;
    % Get latest release tag
    release = git_data.tag_name;
end
