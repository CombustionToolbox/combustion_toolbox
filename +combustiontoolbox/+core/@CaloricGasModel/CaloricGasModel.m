classdef CaloricGasModel < uint8
    % The :mat:class:`CaloricGasModel` class is an enumeration that defines
    % the available caloric gas models for thermodynamic property calculations.
    %
    % The available caloric gas models are:
    %     * `perfect`: Calorically perfect gas model
    %     * `thermallyPerfect`: Thermally perfect gas model
    %     * `imperfect`: Calorically imperfect gas model
    %
    % Example:
    %     model = CaloricGasModel.perfect;

    enumeration
        perfect          (0)
        thermallyPerfect (1)
        imperfect        (2)
    end

    methods

        function obj = setPerfect(obj)
            % Set the caloric gas model to perfect.
            %
            % Example:
            %     model = model.setPerfect()

            obj = combustiontoolbox.core.CaloricGasModel.perfect;
        end

        function obj = setThermallyPerfect(obj)
            % Set the caloric gas model to thermally perfect.
            %
            % Example:
            %     model = model.setThermallyPerfect()

            obj = combustiontoolbox.core.CaloricGasModel.thermallyPerfect;
        end

        function obj = setImperfect(obj)
            % Set the caloric gas model to imperfect.
            %
            % Example:
            %     model = model.setImperfect()

            obj = combustiontoolbox.core.CaloricGasModel.imperfect;
        end

        function value = isPerfect(obj)
            % Check if the caloric gas model is perfect.
            %
            % Returns:
            %     value (logical): True if the model is perfect, false otherwise.
            %
            % Example:
            %     isPerfect = model.isPerfect()

            value = (obj == combustiontoolbox.core.CaloricGasModel.perfect);
        end

        function value = isThermallyPerfect(obj)
            % Check if the caloric gas model is thermally perfect.
            %
            % Returns:
            %     value (logical): True if the model is thermally perfect, false otherwise.
            %
            % Example:
            %     isThermallyPerfect = model.isThermallyPerfect()

            value = (obj == combustiontoolbox.core.CaloricGasModel.thermallyPerfect);
        end

        function value = isImperfect(obj)
            % Check if the caloric gas model is imperfect.
            %
            % Returns:
            %     value (logical): True if the model is imperfect, false otherwise.
            %
            % Example:
            %     isImperfect = model.isImperfect()

            value = (obj == combustiontoolbox.core.CaloricGasModel.imperfect);
        end

    end

    methods (Static)

        function model = fromFlag(FLAG_TCHEM_FROZEN, FLAG_FROZEN)
            % Convert legacy solver flags to CaloricGasModel enumeration.
            %
            % Args:
            %     FLAG_TCHEM_FROZEN (logical): True if calorically perfect gas
            %     FLAG_FROZEN (logical): True if thermally perfect gas
            %
            % Returns:
            %     model (CaloricGasModel): Corresponding caloric gas model
            %
            % Legacy compatibility:
            %     * FLAG_TCHEM_FROZEN = true,  FLAG_FROZEN = false  --- perfect
            %     * FLAG_TCHEM_FROZEN = false, FLAG_FROZEN = true   --- thermallyPerfect
            %     * FLAG_TCHEM_FROZEN = false, FLAG_FROZEN = false  --- imperfect
            %
            % Example:
            %     model = CaloricGasModel.fromFlag(true, false)

            % Check inputs
            if ~islogical(FLAG_TCHEM_FROZEN) || ~islogical(FLAG_FROZEN)
                error('Input flags must be logical values.');
            end

            % Map legacy flags to caloric model
            if FLAG_TCHEM_FROZEN && FLAG_FROZEN
                error('Incompatible flags: Both FLAG_TCHEM_FROZEN and FLAG_FROZEN cannot be true simultaneously.');
            elseif FLAG_TCHEM_FROZEN
                model = combustiontoolbox.core.CaloricGasModel.perfect;
            elseif FLAG_FROZEN
                model = combustiontoolbox.core.CaloricGasModel.thermallyPerfect;
            else
                model = combustiontoolbox.core.CaloricGasModel.imperfect;
            end

        end

    end

end