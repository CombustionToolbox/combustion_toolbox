function [mix3, varargout] = rocket_parameters(mix2, mix3, gravity, varargin)
    % Compute Rocket performance parameters at the throat
    %
    % This method is based on Gordon, S., & McBride, B. J. (1994). NASA reference publication,
    % 1311.
    %
    % Args:
    %     mix2 (struct): Properties of the mixture at the outlet of the chamber
    %     mix3 (struct): Properties of the mixture at the throat
    %     gravity (float): Gravitational acceleration [m/s2]
    %
    % Returns:
    %     mix3 (struct): Properties of the mixture at the throat

    mix3.area_per_mass_flow_rate = area_per_mass_flow_rate(mix3);
    mix3.cstar = characteristic_velocity(mix2, mix3);
    mix3.cf = velocity_relative(mix3) / mix3.cstar;
    [mix3.I_sp, mix3.I_vac] = specific_impulse(mix3, gravity);

    for i = nargin - 3:-1:1

        if ~isempty(varargin{i})
            varargin{i}.cstar = mix3.cstar;
            varargin{i}.cf = velocity_relative(varargin{i}) / varargin{i}.cstar;
            [varargin{i}.I_sp, varargin{i}.I_vac] = specific_impulse(varargin{i}, gravity);
            varargout{i} = varargin{i};
        else
            varargout{i} = [];
        end

    end

end

% SUB-PASS FUNCTIONS
function value = characteristic_velocity(mix2, mix3)
    value = pressure(mix2) * area_per_mass_flow_rate(mix3) * 1e5;
end

function value = area_per_mass_flow_rate(mix)
    value = 1 / (density(mix) * velocity_relative(mix));
end

function [I_sp, I_vac] = specific_impulse(mix, gravity)
    I_sp = velocity_relative(mix) / gravity;
    I_vac = I_sp + (pressure(mix) * area_per_mass_flow_rate(mix) * 1e5) / gravity;
end
