function opts = hjb_options(opts,method)
%HJ B_OPTIONS Common options (no reliance on global OPTSET state).
if nargin<1 || isempty(opts), opts=struct(); end
if nargin<2, method='fd'; end
base=struct('maxit',200,'tol',1e-7,'res_tol',1e-6,'Delta',1000,...
 'showiters',true,'check_generator',true,'generator_tol',1e-9,...
 'infeasible_tol',1e-11,'nres',2,'jump_outside','error',...
 'max_jump_dim',5,'collocation_damping',1,'stop_on_residual',true,...
 'obstacle_tol',1e-6,'lcp_backend','compecon','lcp_tol',1e-7,...
 'qvi_damping',1);
keys=fieldnames(base);
for k=1:numel(keys)
    if ~isfield(opts,keys{k}) || isempty(opts.(keys{k}))
        opts.(keys{k})=base.(keys{k});
    end
end
if opts.maxit<1 || opts.maxit~=floor(opts.maxit) || opts.Delta<=0 || ...
        opts.tol<=0 || opts.res_tol<=0 || opts.collocation_damping<=0 || ...
        opts.collocation_damping>1 || opts.qvi_damping<=0 || opts.qvi_damping>1
    error('hjb:options','Invalid numerical options');
end
if ~ismember(lower(method),{'fd','collocation'})
    error('hjb:method','Unknown solver method');
end
end
