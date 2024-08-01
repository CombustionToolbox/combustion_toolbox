function website_CT()
    % Open Combustion Toolbox's website in default web browser
    
    % Import packages
    import combustiontoolbox.common.Constants
    import combustiontoolbox.utils.openWebsite

    % Definitions
    url = Constants.url.website;
    
    % Open website
    openWebsite(url)
end