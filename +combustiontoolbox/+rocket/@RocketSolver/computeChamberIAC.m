function mix = computeChamberIAC(obj, mix)
    % Compute chemical equilibria at the exit of the chamber (HP) using
    % the Infinite-Area-Chamber (IAC) model
    %
    % This method is based on the method outlined in Gordon, S., & McBride,
    % B. J. (1994). NASA reference publication, 1311.
    %
    % Args:
    %     obj (RocketSolver): RocketSolver object
    %     mix (Mixture): Properties of the initial mixture
    %
    % Returns:
    %     mix (Mixture): Properties of the mixture at the outlet of the chamber
    %
    % Example:
    %     mix = computeChamberIAC(obj, mix1, mix2)

    % Definitions
    obj.equilibriumSolver.problemType = 'HP';
    obj.equilibriumSolver.FLAG_FROZEN = false;

    % Compute chemical equilibria at the exit of the chamber (HP)
    solve(obj.equilibriumSolver, mix);
    
    % Set A_chamber/A_throat
    mix.areaRatio = Inf;
    
    % Restore problemType
    mix.problemType = 'ROCKET_IAC';
end