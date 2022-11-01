function license_content = GPL()
    % Return Combustion Toolbox license
    %
    % Args:
    %     none (empty)
    %
    % Returns:
    %     license_content (str): license content

    license_content = fileread('LICENSE.md');
end