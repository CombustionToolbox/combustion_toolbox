function mix3 = compute_rocket_parameters(mix2, mix3)
    % Compute Rocket performance parameters
    
    mix3.area_per_mass_flow_rate = area_per_mass_flow_rate(mix3);
    mix3.cstar = characteristic_velocity(mix2, mix3);
    mix3.cf = velocity_relative(mix3) / mix3.cstar;
    [mix3.I_sp, mix3.I_vac] = specific_impulse(mix3);
end

% SUB-PASS FUNCTIONS
function value = characteristic_velocity(mix2, mix3)
    value = pressure(mix2) * area_per_mass_flow_rate(mix3) * 1e5;
end

function value = area_per_mass_flow_rate(mix)
    value = 1 / (density(mix) * velocity_relative(mix));
end

function [I_sp, I_vac] = specific_impulse(mix)
    I_sp = velocity_relative(mix);
    I_vac = I_sp + pressure(mix) * area_per_mass_flow_rate(mix) * 1e5;
end