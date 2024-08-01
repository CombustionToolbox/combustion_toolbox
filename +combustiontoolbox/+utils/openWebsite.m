function openWebsite(url)
    % Open website in default web browser
    %   
    %  Example:
    %    openWebsite('https://combustion-toolbox-website.readthedocs.io/')
    
    % Import packages
    import combustiontoolbox.utils.getOS
    
    % Get operating system
    os = getOS();
    
    % Open website
    switch lower(os)
        case 'windows'
            system( sprintf('start %s', url) );
        case 'macos'
            system( sprintf('open %s', url) );
        case 'linux/unix'
            system( sprintf('xdg-open %s', url) );
        otherwise
            error('Operating system not recognized');
    end

end