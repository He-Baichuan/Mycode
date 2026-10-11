function run_setup()
%RUN_SETUP Set up ContinuousTimeDP v4 + check CompEcon dependencies.
root=fileparts(mfilename('fullpath'));
addpath(root,fullfile(root,'examples'),fullfile(root,'tests'));
if exist('fundefn','file')==0 || exist('funbas','file')==0 || ...
        exist('gridmake','file')==0
    error('CompEcon CEtools not on path. Add its location before run_setup');
end
fprintf('ContinuousTimeDP v4 loaded. CompEcon basis routines available.\n');
if exist('lcpsolve','file')~=2
    fprintf('Warning: CompEcon lcpsolve not found: VI/QVI examples require it.\n');
end
fprintf('  HJBsolve:       %s\n',which('HJBsolve'));
fprintf('  scsolve_jump:   %s\n',which('scsolve_jump'));
end
