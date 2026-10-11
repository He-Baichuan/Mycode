function [G,meta]=hjb_intervention_operator(model,grid,W,opts)
%HJB_INTERVENTION_OPERATOR Unified exogenous obstacle and endogenous impulse.
% Model interfaces:
%   model.stop = @(S,j,params{:}) numeric N-by-1, scalar, or numeric N-by-R.
%   model.intervention = @(S,j,W,evalfun,params{:}) struct array with:
%        target: N-by-d points; to: destination regime (scalar or N-vector);
%        reward: scalar or N-vector; feasible (optional): N-vector boolean.
%   Optional direct 'value' candidate is supported for user-optimized targets.
%   evalfun(T,to) evaluates V_to at arbitrary target coordinates.
% If multiple interventions are returned, take their maximum rowwise.
N=grid.N; R=size(W,2); d=grid.d;
G=-inf(N,R); target=nan(N,d,R); toreg=zeros(N,R);
kind=zeros(N,R); choice=zeros(N,R); policy=cell(1,R);
params={}; if isfield(model,'params'), params=model.params; end
outside='error'; if nargin>3 && isfield(opts,'jump_outside'), outside=opts.jump_outside; end
for reg=1:R
    stopvalue=[];
    if isfield(model,'stop') && ~isempty(model.stop)
        stopvalue=model.stop;
    elseif isfield(model,'obstacle') && ~isempty(model.obstacle)
        stopvalue=model.obstacle;
    end
    if ~isempty(stopvalue)
        if isa(stopvalue,'function_handle')
            st=stopvalue(grid.S,reg,params{:});
        else
            st=stopvalue;
            if isvector(st) && numel(st)==R && R>1
                st=st(reg);
            elseif isequal(size(st),[N,R])
                st=st(:,reg);
            end
        end
        if isscalar(st), st=repmat(st,N,1); else, st=st(:); end
        if numel(st)~=N || any(isnan(st))
            error('hjb:obstacle','Stop payoff must be scalar or N-by-1');
        end
        G(:,reg)=st;
        kind(isfinite(st),reg)=1;
    end
    if isfield(model,'intervention') && ~isempty(model.intervention)
        evalfun=@(T,j)hjb_value_eval(grid,W,T,j,outside);
        options=model.intervention(grid.S,reg,W,evalfun,params{:});
        if ~isempty(options) && ~isstruct(options)
            error('hjb:qvi','Intervention callback must return struct array');
        end
        for k=1:numel(options)
            opt=options(k);
            if isfield(opt,'to') && ~isempty(opt.to), j=opt.to;
            else, j=reg; end
            if isscalar(j), j=repmat(j,N,1); else, j=j(:); end
            if numel(j)~=N || any(j<1|j>R|j~=round(j))
                error('hjb:qvi','Invalid intervention destination regime');
            end
            if isfield(opt,'target') && ~isempty(opt.target)
                t=opt.target;
                if ~isequal(size(t),[N,d]), error('hjb:qvi','Intervention target must be N-by-d'); end
            else
                t=grid.S;
            end
            if isfield(opt,'value') && ~isempty(opt.value)
                candidate=local_column(opt.value,N,'value');
            else
                candidate=evalfun(t,j);
                if isfield(opt,'reward') && ~isempty(opt.reward)
                    candidate=candidate+local_column(opt.reward,N,'reward');
                end
            end
            feasible=true(N,1);
            if isfield(opt,'feasible') && ~isempty(opt.feasible)
                feasible=logical(local_column(opt.feasible,N,'feasible'));
            end
            candidate(~feasible)=-inf;
            better=candidate>G(:,reg);
            G(better,reg)=candidate(better);
            kind(better,reg)=2;
            choice(better,reg)=k;
            target(better,:,reg)=t(better,:);
            toreg(better,reg)=j(better);
            if isfield(opt,'action') && ~isempty(opt.action)
                action=opt.action;
                if isscalar(action), action=repmat(action,N,1); end
                if size(action,1)~=N, error('hjb:qvi','Intervention action must have N rows'); end
                if isempty(policy{reg}), policy{reg}=nan(N,size(action,2)); end
                if size(policy{reg},2)~=size(action,2)
                    error('hjb:qvi','All interventions need same action width');
                end
                policy{reg}(better,:)=action(better,:);
            end
        end
    end
end
meta=struct('target',target,'to',toreg,'kind',kind,'choice',choice,'action',{policy});
end
function y=local_column(x,N,name)
if isscalar(x), y=repmat(x,N,1); else, y=x(:); end
if numel(y)~=N || any(isnan(y)) || any(~isfinite(y))
    error('hjb:qvi','%s must be finite scalar or N-vector',name);
end
end
