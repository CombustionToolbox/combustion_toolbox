function header = prettyHeader(txt, varargin)
    % Display a formatted header in the command window.
    %
    % Args:
    %     txt (char): Text to be displayed
    %
    % Optional Args:
    %     % style (char): Style of the header ('unicode', 'ascii', 'banner'). Default is 'unicode'.
    %
    % Returns:
    %     header (char): Formatted header string
    %
    % Example:
    %     prettyHeader('My Header', 'style', 'ascii')
    
    % Default values
    defaultStyle = 'unicode';

    % Parse optional inputs
    p = inputParser;
    addRequired(p, 'txt', @(x) ischar(x) || isstring(x));
    addOptional(p, 'style', defaultStyle, @(x) ischar(x) || isstring(x));
    parse(p, txt, varargin{:});

    % Unpack inputs
    txt = p.Results.txt;
    style = p.Results.style;

    % Create header based on style
    switch lower(style)
        case 'unicode'
            % Explicit mapping for bold alphanumeric characters
            normal = ['A':'Z' 'a':'z' '0':'9' ' '];

            bold = { ...
                '𝗔','𝗕','𝗖','𝗗','𝗘','𝗙','𝗚','𝗛','𝗜','𝗝','𝗞','𝗟','𝗠','𝗡','𝗢','𝗣','𝗤','𝗥','𝗦','𝗧','𝗨','𝗩','𝗪','𝗫','𝗬','𝗭', ...
                '𝗮','𝗯','𝗰','𝗱','𝗲','𝗳','𝗴','𝗵','𝗶','𝗷','𝗸','𝗹','𝗺','𝗻','𝗼','𝗽','𝗾','𝗿','𝘀','𝘁','𝘂','𝘃','𝘄','𝘅','𝘆','𝘇', ...
                '𝟬','𝟭','𝟮','𝟯','𝟰','𝟱','𝟲','𝟳','𝟴','𝟵',' ' ...
            };

            if numel(normal) ~= numel(bold)
                error('Mapping mismatch: normal=%d, bold=%d', numel(normal), numel(bold));
            end

            % Build dictionary
            keys = cellstr(num2cell(normal));
            boldMap = containers.Map(keys, bold);

            % Transform
            header = '';
            for c = char(txt)
                key = string(c);
                if isKey(boldMap, key)
                    header = [header boldMap(key)]; %#ok<AGROW>
                else
                    header = [header c]; %#ok<AGROW>
                end
            end

        case 'ascii'
            header = upper(txt);
            underline = repmat('-', 1, numel(header));
            header = sprintf('%s\n%s', header, underline);

        case 'banner'
            header = sprintf('********** %s **********', upper(txt));

        otherwise
            error('Unknown style: %s. Use "unicode", "ascii", or "banner".', style);
    end

end