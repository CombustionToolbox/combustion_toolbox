function os = getOS()
    % This function returns a char indicating the operating system.
    % 
    % Returns:
    %   os (char): Operating system
    %
    % Example:
    %   os = getOS();
    
    if ispc
        os = 'Windows';
    elseif ismac
        os = 'macOS';
    elseif isunix
        os = 'Linux/Unix';
    else
        os = 'Unknown';
    end
    
end