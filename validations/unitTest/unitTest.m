classdef unitTest < matlab.unittest.TestCase
    % Class that defines a unit test for the Combustion Toolbox.
    % 
    % The unit test compares the results from the Combustion Toolbox with
    % the results from NASA's Chemical with Applications (CEA) code to
    % verify the accuracy of the calculations.
    %
    % Examples:
    %     * results = run(unitTest);
    %     * results = unitTest().run;

    properties (Constant)
        max_rel_error = 0.02; % Max relative error 2%
    end
    
    properties (Access = public)
        database; % Database with PCHIP polynomial fits
    end

    properties (TestParameter)
        temperature = num2cell(200:100:5000);
        vSpecific = num2cell(round([0.2202, 0.4404, 0.6605, 0.8807], 4));
        equivalenceRatio = num2cell(round(0.5:0.1:4, 4));
        velocityPreshock = {505, 3140, 7830, 11862};
        velocityPreshockReflected = {505, 3140, 5169, 5617, 6103, 6632, 7206, 7830};
    end

    methods (Test)
        function obj = unitTest(varargin)
            % Initialization
            
            % Change current directory
            addpath('../../');

            % Unpack inputs
            if nargin
                obj.database = varargin{1}.database;
                return
            end

            % Get database
            obj.database = combustiontoolbox.databases.NasaDatabase();
        end

        function test_TP_CEA_1(testCase, temperature)
            [act_max_rel_error_moles, act_max_rel_error_prop] = run_test_TP_CEA_1(temperature, testCase.database);
            verifyLessThanOrEqual(testCase, act_max_rel_error_moles, testCase.max_rel_error);
            verifyLessThanOrEqual(testCase, act_max_rel_error_prop, testCase.max_rel_error);
        end

        function test_TV_CEA_1(testCase, equivalenceRatio)
            [act_max_rel_error_moles, act_max_rel_error_prop] = run_test_TV_CEA_1(equivalenceRatio, testCase.database);
            verifyLessThanOrEqual(testCase, act_max_rel_error_moles, testCase.max_rel_error);
            verifyLessThanOrEqual(testCase, act_max_rel_error_prop, testCase.max_rel_error);
        end

        function test_TV_CEA_2(testCase, equivalenceRatio)
            [act_max_rel_error_moles, act_max_rel_error_prop] = run_test_TV_CEA_2(equivalenceRatio, testCase.database);
            verifyLessThanOrEqual(testCase, act_max_rel_error_moles, testCase.max_rel_error);
            verifyLessThanOrEqual(testCase, act_max_rel_error_prop, testCase.max_rel_error);
        end

        function test_HP_CEA_1(testCase, equivalenceRatio)
            [act_max_rel_error_moles, act_max_rel_error_prop] = run_test_HP_CEA_1(equivalenceRatio, testCase.database);
            verifyLessThanOrEqual(testCase, act_max_rel_error_moles, testCase.max_rel_error);
            verifyLessThanOrEqual(testCase, act_max_rel_error_prop, testCase.max_rel_error);
        end

        function test_HP_CEA_2(testCase, equivalenceRatio)
            [act_max_rel_error_moles, act_max_rel_error_prop] = run_test_HP_CEA_2(equivalenceRatio, testCase.database);
            verifyLessThanOrEqual(testCase, act_max_rel_error_moles, testCase.max_rel_error);
            verifyLessThanOrEqual(testCase, act_max_rel_error_prop, testCase.max_rel_error);
        end

        function test_EV_CEA_1(testCase, equivalenceRatio)
            [act_max_rel_error_moles, act_max_rel_error_prop] = run_test_EV_CEA_1(equivalenceRatio, testCase.database);
            verifyLessThanOrEqual(testCase, act_max_rel_error_moles, testCase.max_rel_error);
            verifyLessThanOrEqual(testCase, act_max_rel_error_prop, testCase.max_rel_error);
        end
        
        % function test_SV_CEA_1(testCase, vSpecific)
        %     [act_max_rel_error_moles, act_max_rel_error_prop] = run_test_SV_CEA_1(vSpecific, testCase.database);
        %     verifyLessThanOrEqual(testCase, act_max_rel_error_moles, testCase.max_rel_error);
        %     verifyLessThanOrEqual(testCase, act_max_rel_error_prop, testCase.max_rel_error);
        % end
        
        function test_SHOCK_IONIZATION_CEA_1(testCase, velocityPreshock)
            [act_max_rel_error_prop_mix1, act_max_rel_error_prop_mix2] = run_test_SHOCK_IONIZATION_CEA_1(velocityPreshock, testCase.database);
            verifyLessThanOrEqual(testCase, act_max_rel_error_prop_mix1, testCase.max_rel_error);
            verifyLessThanOrEqual(testCase, act_max_rel_error_prop_mix2, testCase.max_rel_error);
        end
        
        function test_SHOCK_R_IONIZATION_CEA_1(testCase, velocityPreshockReflected)
            [act_max_rel_error_prop_mix1, act_max_rel_error_prop_mix2] = run_test_SHOCK_R_IONIZATION_CEA_1(velocityPreshockReflected, testCase.database);
            verifyLessThanOrEqual(testCase, act_max_rel_error_prop_mix1, testCase.max_rel_error);
            verifyLessThanOrEqual(testCase, act_max_rel_error_prop_mix2, testCase.max_rel_error);
        end
        
        function test_DET_CEA_1(testCase, equivalenceRatio)
            [act_max_rel_error_moles, act_max_rel_error_prop] = run_test_DET_CEA_1(equivalenceRatio, testCase.database);
            verifyLessThanOrEqual(testCase, act_max_rel_error_moles, testCase.max_rel_error);
            verifyLessThanOrEqual(testCase, act_max_rel_error_prop, testCase.max_rel_error);
        end

        function test_DET_CEA_2(testCase, equivalenceRatio)

            % CEA does not converge for phi = 3.9
            if equivalenceRatio == 3.9
                equivalenceRatio = 4;
                return
            end

            [act_max_rel_error_moles, act_max_rel_error_prop] = run_test_DET_CEA_2(equivalenceRatio, testCase.database);
            verifyLessThanOrEqual(testCase, act_max_rel_error_moles, testCase.max_rel_error);
            verifyLessThanOrEqual(testCase, act_max_rel_error_prop, testCase.max_rel_error);
        end
         
        function test_ROCKET_CEA_1(testCase, equivalenceRatio)

            % CEA does not converge for phi = 4
            if equivalenceRatio > 3.9
                equivalenceRatio = 3.9;
            end

            [act_max_rel_error_moles, act_max_rel_error_prop] = run_test_ROCKET_CEA_1(equivalenceRatio, testCase.database);
            verifyLessThanOrEqual(testCase, act_max_rel_error_moles, testCase.max_rel_error);
            verifyLessThanOrEqual(testCase, act_max_rel_error_prop, testCase.max_rel_error);
        end

    end

end
