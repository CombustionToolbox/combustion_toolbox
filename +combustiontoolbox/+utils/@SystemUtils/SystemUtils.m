classdef SystemUtils < handle
    % SystemUtils - A class containing utility functions for system operations
    %
    % Examples:
    %   * SystemUtils.openWebsite('https://combustion-toolbox-website.readthedocs.io/');
    %   * SystemUtils.openWebsite(SystemUtils.url.documentation);
    %   * os = SystemUtils.getOS();

    properties (Constant)
        url = getUrl() % URLs
    end

    methods (Static)
        function os = getOS()
            % This function returns a char indicating the operating system.
            % 
            % Returns:
            %   os (char): Operating system
            %
            % Example:
            %   os = SystemUtils.getOS();
            
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

        function openWebsite(url)
            % Open website in default web browser
            %   
            %  Example:
            %    SystemUtils.openWebsite('https://combustion-toolbox-website.readthedocs.io/')
            
            % Import packages
            import combustiontoolbox.utils.SystemUtils
            
            % Get operating system
            os = SystemUtils.getOS();
            
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

        function websiteCT()
            % Open Combustion Toolbox's website in default web browser
            
            % Import packages
            import combustiontoolbox.utils.SystemUtils
        
            % Definitions
            url = SystemUtils.url.website;
            
            % Open website
            SystemUtils.openWebsite(url)
        end

    end

end

% SUB-PASS FUNCTIONS
function url = getUrl()
    url.repository = 'https://github.com/CombustionToolbox/combustion_toolbox';
    url.website = 'https://combustion-toolbox-website.readthedocs.io/';
    url.documentation = 'https://combustion-toolbox-website.readthedocs.io/en/latest/documentation/api/index.html';
    url.publications = 'https://combustion-toolbox-website.readthedocs.io/en/latest/publications.html';
    url.websiteUC3M = 'http://fluidosuc3m.es/';
    url.websiteACuadra = 'https://acuadralara.com/';
    url.mailACuadra = 'mailto:acuadra@ing.uc3m.es';
    url.githubACuadra = 'https://github.com/AlbertoCuadra';
    url.researchgateACuadra = 'https://www.researchgate.net/profile/Alberto_Cuadra_Lara';
    url.linkedinACuadra = 'https://www.linkedin.com/in/albertocuadralara/';
    url.websiteMVera = 'http://fluidosuc3m.es/people/mvcoello/';
    url.researchgateMVera = 'https://www.researchgate.net/profile/Marcos_Vera';
    url.linkedinMVera = 'https://www.linkedin.com/in/marcos-vera-coello-05b3a643/?originalSubdomain=es';
    url.websiteCHuete = 'http://fluidosuc3m.es/people/chuete/';
    url.researchgateCHuete = 'https://www.researchgate.net/profile/Cesar-Huete';
    url.linkedinCHuete = 'https://www.linkedin.com/in/cesarhuete/';
    url.contributors = 'https://github.com/CombustionToolbox/combustion_toolbox/blob/master/CONTRIBUTORS.md';
end