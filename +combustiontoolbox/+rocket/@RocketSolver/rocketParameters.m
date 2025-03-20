function [mix3, varargout] = rocketParameters(mix2, mix3, varargin)
    % Compute Rocket performance parameters at the throat
    %
    % This method is based on the method outlined in Gordon, S., & McBride,
    % B. J. (1994). NASA reference publication, 1311.
    %
    % Args:
    %     mix2 (Mixture): Properties of the mixture at the outlet of the chamber
    %     mix3 (Mixture): Properties of the mixture at the throat
    %
    % Optional Args:
    %     * mixi (Mixture): Properties of the mixture at the the exit (obj.FLAG_SUBSONIC = false) or at the combustion chamber (obj.FLAG_SUBSONIC = true)
    %
    % Returns:
    %     mix3 (Mixture): Properties of the mixture at the throat

    % Definitions
    gravity = combustiontoolbox.common.Constants.G;

    % Compute rocket parameters
    mix3.cstar = characteristicVelocity(mix2, mix3);
    mix3.cf = mix3.u / mix3.cstar;
    [mix3.I_sp, mix3.I_vac] = specificImpulse(mix3, gravity);

    for i = nargin-2:-1:1
        
        if ~isempty(varargin{i})
            varargin{i}.cstar = mix3.cstar;
            varargin{i}.cf = varargin{i}.u / varargin{i}.cstar;
            [varargin{i}.I_sp, varargin{i}.I_vac] = specificImpulse(varargin{i}, gravity);
        end

        varargout{i} = varargin{i};
    end

end

% SUB-PASS FUNCTIONS
function value = characteristicVelocity(mix2, mix3)
    % Compute the characteristic velocity
    value = mix2.p * areaPerMassFlowRate(mix3) * 1e5;
end

function value = areaPerMassFlowRate(mix)
    % Compute the area per unit mass flow rate
    value = 1 / (mix.rho * mix.u);
end

function [I_sp, I_vac] = specificImpulse(mix, gravity)
    % Compute specific impulse values (sea level and vacuum)
    I_sp = mix.u / gravity;
    I_vac = I_sp + (mix.p * areaPerMassFlowRate(mix) * 1e5) / gravity;
end
