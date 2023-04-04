classdef unit_test < matlab.unittest.TestCase
    % Class that defines a unit test for the Combustion Toolbox.
    % 
    % The unit test compares the results from the Combustion Toolbox with
    % the results from NASA's Chemical with Applications (CEA) code to
    % verify the accuracy of the calculations.
    %
    % Local test:
    %     * results = run(unit_test);
    %     * results = unit_test().run;

    properties (Constant)
        max_rel_error = 0.02; % Max relative error 2%
    end
    
    properties (Access = public)
        DB; % Database with PCHIP polynomial fits
        DB_master; % Database    
    end

    properties (TestParameter)
        equivalence_ratio = {0.5, 1, 2, 4};
        velocity_preshock = {505, 3140, 7830, 11862};
        velocity_preshock_R = {505, 3140, 5169, 5617, 6103, 6632, 7206, 7830};
    end

    methods (Test)
        function obj = unit_test(varargin)
            % Initialization

            % Unpack inputs
            if nargin
                obj.DB = varargin{1}.DB;
                obj.DB_master = varargin{1}.DB_master;
                return
            end

            % Get databases
            self = App();
            obj.DB = self.DB;
            obj.DB_master = self.DB;
        end

        function test_HP_CEA_1(testCase, equivalence_ratio)
            [act_max_rel_error_moles, act_max_rel_error_prop] = run_test_HP_CEA_1(equivalence_ratio, testCase.DB, testCase.DB_master);
            verifyLessThanOrEqual(testCase, act_max_rel_error_moles, testCase.max_rel_error);
            verifyLessThanOrEqual(testCase, act_max_rel_error_prop, testCase.max_rel_error);
        end

        function test_HP_CEA_2(testCase, equivalence_ratio)
            [act_max_rel_error_moles, act_max_rel_error_prop] = run_test_HP_CEA_2(equivalence_ratio, testCase.DB, testCase.DB_master);
            verifyLessThanOrEqual(testCase, act_max_rel_error_moles, testCase.max_rel_error);
            verifyLessThanOrEqual(testCase, act_max_rel_error_prop, testCase.max_rel_error);
        end

        function test_EV_CEA_1(testCase, equivalence_ratio)
            [act_max_rel_error_moles, act_max_rel_error_prop] = run_test_EV_CEA_1(equivalence_ratio, testCase.DB, testCase.DB_master);
            verifyLessThanOrEqual(testCase, act_max_rel_error_moles, testCase.max_rel_error);
            verifyLessThanOrEqual(testCase, act_max_rel_error_prop, testCase.max_rel_error);
        end

        function test_SHOCK_IONIZATION_CEA_1(testCase, velocity_preshock)
            [act_max_rel_error_prop_mix1, act_max_rel_error_prop_mix2] = run_test_SHOCK_IONIZATION_CEA_1(velocity_preshock, testCase.DB, testCase.DB_master);
            verifyLessThanOrEqual(testCase, act_max_rel_error_prop_mix1, testCase.max_rel_error);
            verifyLessThanOrEqual(testCase, act_max_rel_error_prop_mix2, testCase.max_rel_error);
        end

        function test_SHOCK_R_IONIZATION_CEA_1(testCase, velocity_preshock_R)
            [act_max_rel_error_prop_mix1, act_max_rel_error_prop_mix2] = run_test_SHOCK_R_IONIZATION_CEA_1(velocity_preshock_R, testCase.DB, testCase.DB_master);
            verifyLessThanOrEqual(testCase, act_max_rel_error_prop_mix1, testCase.max_rel_error);
            verifyLessThanOrEqual(testCase, act_max_rel_error_prop_mix2, testCase.max_rel_error);
        end

        function test_DET_CEA_1(testCase, equivalence_ratio)
            [act_max_rel_error_moles, act_max_rel_error_prop] = run_test_DET_CEA_1(equivalence_ratio, testCase.DB, testCase.DB_master);
            verifyLessThanOrEqual(testCase, act_max_rel_error_moles, testCase.max_rel_error);
            verifyLessThanOrEqual(testCase, act_max_rel_error_prop, testCase.max_rel_error);
        end

        % function test_DET_CEA_2(testCase, equivalence_ratio)
        %     [act_max_rel_error_moles, act_max_rel_error_prop] = run_test_DET_CEA_2(equivalence_ratio);
        %     verifyLessThanOrEqual(testCase, act_max_rel_error_moles, testCase.max_rel_error);
        %     verifyLessThanOrEqual(testCase, act_max_rel_error_prop, testCase.max_rel_error);
        % end

    end

end
