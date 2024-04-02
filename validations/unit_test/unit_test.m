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
    end

    properties (TestParameter)
        temperature = {200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800};
        density = {4.5418, 2.2709, 1.5139, 1.1355};
        equivalence_ratio = {0.5, 1, 2, 4};
        velocity_preshock = {505, 3140, 7830, 11862};
        velocity_preshock_R = {505, 3140, 5169, 5617, 6103, 6632, 7206, 7830};
    end

    methods (Test)
        function obj = unit_test(varargin)
            % Initialization
            
            % Change current directory
            addpath('../../');

            % Unpack inputs
            if nargin
                obj.DB = varargin{1}.DB;
                return
            end

            % Get database
            obj.DB = combustiontoolbox.databases.NasaDatabase();
        end

        function test_TP_CEA_1(testCase, temperature)
            [act_max_rel_error_moles, act_max_rel_error_prop] = run_test_TP_CEA_1(temperature, testCase.DB);
            verifyLessThanOrEqual(testCase, act_max_rel_error_moles, testCase.max_rel_error);
            verifyLessThanOrEqual(testCase, act_max_rel_error_prop, testCase.max_rel_error);
        end

        function test_HP_CEA_1(testCase, equivalence_ratio)
            [act_max_rel_error_moles, act_max_rel_error_prop] = run_test_HP_CEA_1(equivalence_ratio, testCase.DB);
            verifyLessThanOrEqual(testCase, act_max_rel_error_moles, testCase.max_rel_error);
            verifyLessThanOrEqual(testCase, act_max_rel_error_prop, testCase.max_rel_error);
        end

        function test_HP_CEA_2(testCase, equivalence_ratio)
            [act_max_rel_error_moles, act_max_rel_error_prop] = run_test_HP_CEA_2(equivalence_ratio, testCase.DB);
            verifyLessThanOrEqual(testCase, act_max_rel_error_moles, testCase.max_rel_error);
            verifyLessThanOrEqual(testCase, act_max_rel_error_prop, testCase.max_rel_error);
        end

        % function test_EV_CEA_1(testCase, equivalence_ratio)
        %     [act_max_rel_error_moles, act_max_rel_error_prop] = run_test_EV_CEA_1(equivalence_ratio, testCase.DB, testCase.DB_master);
        %     verifyLessThanOrEqual(testCase, act_max_rel_error_moles, testCase.max_rel_error);
        %     verifyLessThanOrEqual(testCase, act_max_rel_error_prop, testCase.max_rel_error);
        % end
        % 
        % function test_SV_CEA_1(testCase, density)
        %     [act_max_rel_error_moles, act_max_rel_error_prop] = run_test_SV_CEA_1(density, testCase.DB, testCase.DB_master);
        %     verifyLessThanOrEqual(testCase, act_max_rel_error_moles, testCase.max_rel_error);
        %     verifyLessThanOrEqual(testCase, act_max_rel_error_prop, testCase.max_rel_error);
        % end
        % 
        % function test_SHOCK_IONIZATION_CEA_1(testCase, velocity_preshock)
        %     [act_max_rel_error_prop_mix1, act_max_rel_error_prop_mix2] = run_test_SHOCK_IONIZATION_CEA_1(velocity_preshock, testCase.DB, testCase.DB_master);
        %     verifyLessThanOrEqual(testCase, act_max_rel_error_prop_mix1, testCase.max_rel_error);
        %     verifyLessThanOrEqual(testCase, act_max_rel_error_prop_mix2, testCase.max_rel_error);
        % end
        % 
        % function test_SHOCK_R_IONIZATION_CEA_1(testCase, velocity_preshock_R)
        %     [act_max_rel_error_prop_mix1, act_max_rel_error_prop_mix2] = run_test_SHOCK_R_IONIZATION_CEA_1(velocity_preshock_R, testCase.DB, testCase.DB_master);
        %     verifyLessThanOrEqual(testCase, act_max_rel_error_prop_mix1, testCase.max_rel_error);
        %     verifyLessThanOrEqual(testCase, act_max_rel_error_prop_mix2, testCase.max_rel_error);
        % end
        % 
        % function test_DET_CEA_1(testCase, equivalence_ratio)
        %     [act_max_rel_error_moles, act_max_rel_error_prop] = run_test_DET_CEA_1(equivalence_ratio, testCase.DB, testCase.DB_master);
        %     verifyLessThanOrEqual(testCase, act_max_rel_error_moles, testCase.max_rel_error);
        %     verifyLessThanOrEqual(testCase, act_max_rel_error_prop, testCase.max_rel_error);
        % end
        % 
        % % function test_DET_CEA_2(testCase, equivalence_ratio)
        % %     [act_max_rel_error_moles, act_max_rel_error_prop] = run_test_DET_CEA_2(equivalence_ratio);
        % %     verifyLessThanOrEqual(testCase, act_max_rel_error_moles, testCase.max_rel_error);
        % %     verifyLessThanOrEqual(testCase, act_max_rel_error_prop, testCase.max_rel_error);
        % % end
        % 
        % function test_ROCKET_CEA_1(testCase, equivalence_ratio)
        % 
        %     % CEA does not converge for phi = 4
        %     if equivalence_ratio > 3.9
        %         equivalence_ratio = 3.9;
        %     end
        % 
        %     [act_max_rel_error_moles, act_max_rel_error_prop] = run_test_ROCKET_CEA_1(equivalence_ratio, testCase.DB, testCase.DB_master);
        %     verifyLessThanOrEqual(testCase, act_max_rel_error_moles, testCase.max_rel_error);
        %     verifyLessThanOrEqual(testCase, act_max_rel_error_prop, testCase.max_rel_error);
        % end

    end

end
