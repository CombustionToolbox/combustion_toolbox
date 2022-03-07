% Routine that checks a set of test of the code 
 
% Run the tests and fail the build if any of the tests fails 
Combustion_Toolbox_setPath()
testCase = app_test;
results = testCase.run; 
nfailed = nnz([results.Failed]);
assert(nfailed == 0,[num2str(nfailed) ' test(s) failed.'])