classdef Constants < handle
    % Class with constants data
    %
    % Attributes:
    %     R0 (float): Universal gas constant [J/(K mol)]
    %     G (float): Standard gravity [m/s2]
    %     NA (float): Avogadro's number [molecule/mol]
    %     release (char): Release of the Combustion Toolbox
    %     date (char): Date of the release
    %     url (struct): URLs
    %
    % Examples:
    %     * R0 = Constants.R0
    %     * g = Constants.G
    %     * release = Constants.release
    
    properties (Constant)
        R0      = 8.31446261815324 % Universal gas constant [J/(K mol)]
        G       = 9.80665          % Standard gravity [m/s2]
        NA      = 6.0221415e23     % Avogadro's number [molecule/mol]
        release = 'v1.1.0'         % Release version
        date    = '07 Jun 2024'    % Release date
        url     = getUrl()         % URLs
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