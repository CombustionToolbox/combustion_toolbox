function mix = rocketChamberIAC(obj, mix, mix_guess)
    % Compute chemical equilibria at the exit of the chamber (HP) using
    % the Infinite-Area-Chamber (IAC) model
    %
    % This method is based on the method outlined in Gordon, S., & McBride,
    % B. J. (1994). NASA reference publication, 1311.
    %
    % Args:
    %     obj (RocketSolver): RocketSolver object
    %     mix (Mixture): Properties of the initial mixture
    %     mix_guess (Mixture): Properties of the mixture at the outlet of the chamber (previous calculation)
    %
    % Returns:
    %     mix (Mixture): Properties of the mixture at the outlet of the chamber
    %
    % Example:
    %     mix = rocketChamberIAC(obj, mix, mix_guess)
    
    % Temporal value
    TEMP_FLAG_FROZEN = obj.equilibriumSolver.FLAG_FROZEN;

    % Definitions
    obj.equilibriumSolver.problemType = 'HP';
    obj.equilibriumSolver.FLAG_FROZEN = false;

    % Compute chemical equilibria at the exit of the chamber (HP)
    if isempty(mix_guess)
        solve(obj.equilibriumSolver, mix);
    else
        solve(obj.equilibriumSolver, mix, mix_guess);
    end

    % Set areaRatio = areaChamber / areaThroat
    mix.areaRatio = Inf;
    
    % Restore problemType and FLAG_FROZEN values;
    mix.problemType = 'ROCKET_IAC';
    obj.equilibriumSolver.FLAG_FROZEN = TEMP_FLAG_FROZEN;
end