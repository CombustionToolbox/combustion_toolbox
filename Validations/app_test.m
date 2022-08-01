classdef app_test < matlab.unittest.TestCase
    % Test solutions from Combustion Toolbox against NASA's Chemical with
    % Applications (CEA) code
    %
    % Local test:
    % testCase = app_test;
    % results = testCase.run;

    properties(Constant)
        max_rel_error = 0.02; % Max relative error 2%
    end

    properties(TestParameter)
        equivalence_ratio = {0.5, 1, 2, 4};
        velocity_preshock = {505, 3140, 7830, 11862};
        velocity_preshock_R = {505, 3140, 5169, 5617, 6103, 6632, 7206, 7830};
    end

    methods (Test)
        function test_HP_CEA_1(testCase, equivalence_ratio)
            [act_max_rel_error_moles, act_max_rel_error_prop] = run_test_HP_CEA_1(equivalence_ratio);
            verifyLessThanOrEqual(testCase, act_max_rel_error_moles, testCase.max_rel_error);
            verifyLessThanOrEqual(testCase, act_max_rel_error_prop, testCase.max_rel_error);
        end
        function test_HP_CEA_2(testCase, equivalence_ratio)
            [act_max_rel_error_moles, act_max_rel_error_prop] = run_test_HP_CEA_2(equivalence_ratio);
            verifyLessThanOrEqual(testCase, act_max_rel_error_moles, testCase.max_rel_error);
            verifyLessThanOrEqual(testCase, act_max_rel_error_prop, testCase.max_rel_error);
        end
        function test_EV_CEA_1(testCase, equivalence_ratio)
            [act_max_rel_error_moles, act_max_rel_error_prop] = run_test_EV_CEA_1(equivalence_ratio);
            verifyLessThanOrEqual(testCase, act_max_rel_error_moles, testCase.max_rel_error);
            verifyLessThanOrEqual(testCase, act_max_rel_error_prop, testCase.max_rel_error);
        end
        function test_SHOCK_IONIZATION_CEA_1(testCase, velocity_preshock)
            [act_max_rel_error_prop_mix1, act_max_rel_error_prop_mix2] = run_test_SHOCK_IONIZATION_CEA_1(velocity_preshock);
            verifyLessThanOrEqual(testCase, act_max_rel_error_prop_mix1, testCase.max_rel_error);
            verifyLessThanOrEqual(testCase, act_max_rel_error_prop_mix2, testCase.max_rel_error);
        end
        function test_SHOCK_R_IONIZATION_CEA_1(testCase, velocity_preshock_R)
            [act_max_rel_error_prop_mix1, act_max_rel_error_prop_mix2] = run_test_SHOCK_R_IONIZATION_CEA_1(velocity_preshock_R);
            verifyLessThanOrEqual(testCase, act_max_rel_error_prop_mix1, testCase.max_rel_error);
            verifyLessThanOrEqual(testCase, act_max_rel_error_prop_mix2, testCase.max_rel_error);
        end
        function test_DET_CEA_1(testCase, equivalence_ratio)
            [act_max_rel_error_moles, act_max_rel_error_prop] = run_test_DET_CEA_1(equivalence_ratio);
            verifyLessThanOrEqual(testCase, act_max_rel_error_moles, testCase.max_rel_error);
            verifyLessThanOrEqual(testCase, act_max_rel_error_prop, testCase.max_rel_error);
        end
%         function test_DET_CEA_2(testCase, equivalence_ratio)
%             [act_max_rel_error_moles, act_max_rel_error_prop] = run_test_DET_CEA_2(equivalence_ratio);
%             verifyLessThanOrEqual(testCase, act_max_rel_error_moles, testCase.max_rel_error);
%             verifyLessThanOrEqual(testCase, act_max_rel_error_prop, testCase.max_rel_error);
%         end
    end
end