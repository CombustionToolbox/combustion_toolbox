fun1 = fullfile('combustion_toolbox','*.m');
fun2 = fullfile('combustion_toolbox','Databases','*.m');
fun3 = fullfile('combustion_toolbox','Display','*.m');
fun4 = fullfile('combustion_toolbox','Export','*.m');
fun5 = fullfile('combustion_toolbox','GUI','*.m');
fun6 = fullfile('combustion_toolbox','NASA_database','*.m');
fun7 = fullfile('combustion_toolbox','Settings','*.m');
fun8 = fullfile('combustion_toolbox','Solver','Chemical Equilibrium','*.m');
fun9 = fullfile('combustion_toolbox','Solver','Functions','*.m');
fun10 = fullfile('combustion_toolbox','Solver','Shocks and detonations','*.m');

pcode(fun1,fun2,fun3,fun4,fun5,fun6,fun7,fun8,fun9,fun10,'-inplace')
