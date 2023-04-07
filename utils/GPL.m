function license_content = GPL()
    % Return Combustion Toolbox license
    %
    % Returns:
    %   license_content (char): The license text
    %
    % Example:
    %   license_content = GPL()

    license_content = fileread('LICENSE.md');
end