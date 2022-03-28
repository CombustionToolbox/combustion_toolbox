function self = set_prop(self, varargin)
    % Assign property values to the respective variables
    %
    % Args:
    %     self (struct): Data of the mixture, conditions, and databases
    %
    % Optional Args:
    %     field (str):   Fieldname
    %     value (float): Value/s to assing in the field
    %
    % Returns:
    %     self (struct): Data of the mixture, conditions, and databases

    n = round(nargin/2 - 1);
    try
        for i=1:n
            self = set_prop_common(self, varargin{2*i - 1}, varargin{2*i});
        end
    catch
        error('Error number inputs @set_prop');
    end
end

% SUB-PASS FUNCTIONS
function self = set_prop_common(self, name, value)
    % Assign property values to the given variable name

    % If the value is a string, convert it to a float - (GUI)
    if value(1) == '['
        value = sscanf(value, '[%f:%f:%f]');
        value = value(1):value(2):value(3);
    elseif ischar(value)
        value = sscanf(value,'%f');
    end
    % Define flags
    if length(value) > 1
        self.Misc.FLAGS_PROP.(name) = true;
    else
        self.Misc.FLAGS_PROP.(name) = false;
    end
    % Set value
    self.PD.(name).value = value;
end