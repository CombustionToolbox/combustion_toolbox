function self = set_prop(varargin)
    self = varargin{1};
    n = round(nargin/2 - 1);
    try
        for i=1:n
            self = set_prop_common(self, varargin{2*i}, varargin{2*i + 1});
        end
    catch
        error('Error number inputs @set_prop');
    end
end

% SUB-PASS FUNCTIONS
function self = set_prop_common(self, name, value)
    if length(value) > 1
        self.Misc.FLAGS_PROP.(name) = true;
    else
        self.Misc.FLAGS_PROP.(name) = false;
    end
    self.PD.(name).value = value;
end