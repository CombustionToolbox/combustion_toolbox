function website_CT()
    % Open Combustion Toolbox's website in default web browser
    
    % Import packages
    import combustiontoolbox.utils.SystemUtils

    % Definitions
    url = SystemUtils.url.website;
    
    % Open website
    SystemUtils.openWebsite(url)
end