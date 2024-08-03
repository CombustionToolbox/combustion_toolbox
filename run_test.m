% Routine that checks a set of unit tests of the code
% 
% Note:
%     Stop the build if any tests failed

% Set the path for the Combustion Toolbox
INSTALL('install', 'path');

% Run the tests and store the results
results = unitTest().run;

% Count the number of failed tests
N_failed = nnz([results.Failed]);

% Stop the build if any tests failed
assert(N_failed == 0, sprintf('%d test(s) failed.', N_failed));