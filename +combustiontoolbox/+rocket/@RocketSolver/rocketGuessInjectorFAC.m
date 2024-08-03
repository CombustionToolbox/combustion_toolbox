function pressure_inf = rocketGuessInjectorFAC(pressure_inj, areaRatioChamber)
    % Compute pressure guess [bar] assuming an Infinite-Area-Chamber (IAC)
    % for the Finite-Area-Chamber model (FAC)
    %
    % Args:
    %     pressure_inj (float): Pressure at the injector [bar]
    %     areaRatioChamber (float): Area chamber / Area throat
    %
    % Returns:
    %     pressure (float): Pressure at the throat [bar]
    %
    % Example:
    %     pressure = rocketGuessInjectorFAC(pressure_inj, areaRatioChamber)

    pressure_inf = pressure_inj * (1.0257 - 1.2318 * areaRatioChamber) / (1 - 1.26505 * areaRatioChamber);
end