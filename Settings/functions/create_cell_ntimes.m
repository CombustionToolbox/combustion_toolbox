function C = create_cell_ntimes(varargin)
    % Create cell array with the same item n-times
    if nargin > 2
        value = varargin{1};
        C = varargin{3};
    elseif nargin > 1
        value = varargin{1};
        n = varargin{2};
        C = cell(1, n);
    else
        error('Error sub-pass fuinction @create_cell_ntimes inside @guiReactantsValueChanged');
    end
    % Set value
    C(:) = {value};
end