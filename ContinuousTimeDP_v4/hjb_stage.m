function out=hjb_stage(model,G,W,opts)
%HJB_STAGE Shared policy-candidate evaluation for FD and collocation.
% model.candidates(S,D,W,reg,params{:}) returns struct array with
% fields x, f, g, sigma(optional), feasible(optional), q(optional), name.
% Scores every candidate by its FULL discretized Hamiltonian, then selects
% rowwise maximizer. Candidate generation must explore global optima.
N=G.N; d=G.d; R=size(W,2);
if ~isfield(model,'candidates') || ~isa(model.candidates,'function_handle')
    error('hjb:callback','model.candidates callback required');
end
params={}; if isfield(model,'params'), params=model.params; end
boundary=repmat({'state_constraint'},1,d);
if isfield(model,'boundary') && ~isempty(model.boundary)
    boundary=model.boundary;
    if ischar(boundary), boundary=repmat({boundary},1,d); end
end
if numel(boundary)~=d, error('hjb:boundary','model.boundary must specify each dimension'); end
if strcmp(G.method,'collocation'), C=G.P0\W; state=C(:); else, C=[]; state=W(:); end
A=sparse(N*R,N*R);
fval=zeros(N,R); pol=cell(1,R); branch=zeros(N,R); rho=zeros(N,R);
candscores=cell(1,R); candidates_count=zeros(1,R);
advpol=cell(1,R); advbranch=zeros(N,R);
for reg=1:R
    rho(:,reg)=hjb_eval_rho(model,G.S,reg);
    if strcmp(G.method,'fd')
        D=hjb_fd_derivatives(G,W(:,reg));
    else
        D.grad=zeros(N,d);
        for k=1:d, D.grad(:,k)=G.P1{k}*C(:,reg); end
        D.forward=D.grad; D.backward=D.grad; D.central=D.grad;
    end
    D.value=W(:,reg); D.all_values=W; D.regime=reg;
    cands=model.candidates(G.S,D,W,reg,params{:});
    if ~isstruct(cands) || isempty(cands)
        error('hjb:policy','Candidate callback returned no controls');
    end
    % For robust problems evaluate nature's alternatives for each controller.
    % Then choose max_{controller} min_{nature} Hamiltonian rowwise.
    ncontroller=numel(cands); candidates_count(reg)=ncontroller;
    expanded={}; group=[]; adversary_x={}; adversary_index=[];
    for ctrl=1:ncontroller
        if isfield(model,'adversary') && ~isempty(model.adversary)
            variants=model.adversary(G.S,D,W,reg,cands(ctrl),params{:});
            if ~isstruct(variants) || isempty(variants)
                error('hjb:adversary','Adversary callback must return struct variants');
            end
        else
            variants=struct('x',[]);
        end
        for nature=1:numel(variants)
            base=cands(ctrl); alt=variants(nature);
            % An adversary's x means its own distortion, not controller x.
            if isfield(alt,'x'), ax=alt.x; else, ax=[]; end
            overrides={'f','g','sigma','q','jump','feasible'};
            for o=1:numel(overrides)
                key=overrides{o};
                if isfield(alt,key) && ~isempty(alt.(key))
                    base.(key)=alt.(key);
                end
            end
            expanded{end+1}=base; %#ok<AGROW>
            adversary_x{end+1}=ax; %#ok<AGROW>
            group(end+1)=ctrl; %#ok<AGROW>
            adversary_index(end+1)=nature; %#ok<AGROW>
        end
    end
    nc=numel(expanded);
    scores=-inf(N,nc); Llist=cell(1,nc); Flist=cell(1,nc); Xlist=cell(1,nc);
    for t=1:nc
        candidate=expanded{t};
        if ~isfield(candidate,'f') || ~isfield(candidate,'g') || ~isfield(candidate,'x')
            error('hjb:policy','Candidate needs x, f, g');
        end
        fi=hjb_column(candidate.f,N,'f');
        mu=candidate.g;
        if d==1 && isvector(mu) && numel(mu)==N, mu=mu(:); end
        if ~isequal(size(mu),[N,d]) || any(~isfinite(mu(:)))
            error('hjb:policy','Candidate drift must be N-by-d');
        end
        xp=candidate.x;
        if isvector(xp) && numel(xp)==N, xp=xp(:); end
        if size(xp,1)~=N || any(~isfinite(xp(:)))
            error('hjb:policy','Candidate controls must be N-by-p finite');
        end
        sig=[];
        if isfield(candidate,'sigma') && ~isempty(candidate.sigma), sig=candidate.sigma; end
        if ~isempty(sig)
            if d==1 && isvector(sig), sig=reshape(sig,N,1,1); end
            if ismatrix(sig) && d>1 && isequal(size(sig),[N,d])
                sig=reshape(sig,N,d,1);
            end
            if size(sig,1)~=N || size(sig,2)~=d || any(~isfinite(sig(:)))
                error('hjb:policy','Candidate sigma must be N-by-d-by-q');
            end
        end
        feasible=true(N,1);
        if isfield(candidate,'feasible') && ~isempty(candidate.feasible)
            feasible=logical(hjb_column(candidate.feasible,N,'feasible'));
        end
        safeMu=mu; safeSig=sig;
        % Do not let rejected (infeasible) candidate rows poison the FD
        % generator; rejecting any row still requires some valid candidate.
        for k=1:d
            if strcmpi(boundary{k},'state_constraint')
                left=G.subs(:,k)==1; right=G.subs(:,k)==G.n(k);
                violation=(left & mu(:,k)<-opts.infeasible_tol) | ...
                          (right & mu(:,k)>opts.infeasible_tol);
                if ~isempty(sig)
                    sig2=sum(reshape(sig(:,k,:),N,[]).^2,2);
                    violation=violation | ((left|right)&sig2>opts.infeasible_tol);
                    safeSig(violation,k,:)=0;
                end
                feasible=feasible & ~violation;
                safeMu(violation,k)=0;
            end
        end
        if strcmp(G.method,'fd')
            Lop=hjb_fd_operator(G,safeMu,safeSig,boundary,opts);
        else
            Lop=hjb_collocation_operator(G,mu,sig);
        end
        AA=sparse(N,N*R);
        jcols=(reg-1)*N+(1:N);
        AA(:,jcols)=AA(:,jcols)+Lop;
        q=hjb_policy_q(model,candidate,G.S,xp,reg,R,params);
        for dest=1:R
            if dest==reg, continue; end
            rate=q(:,dest);
            if ~any(rate), continue; end
            zcols=(dest-1)*N+(1:N);
            AA(:,zcols)=AA(:,zcols)+spdiags(rate,0,N,N)*G.P0;
            AA(:,jcols)=AA(:,jcols)-spdiags(rate,0,N,N)*G.P0;
        end
        channels=[];
        if isfield(candidate,'jump') && ~isempty(candidate.jump)
            channels=candidate.jump;
            if isa(channels,'function_handle'), channels=channels(G.S,xp,reg,params{:}); end
        elseif isfield(model,'jump') && ~isempty(model.jump)
            channels=model.jump(G.S,xp,reg,params{:});
        end
        if ~isempty(channels)
            for ch=1:numel(channels)
                item=channels(ch);
                if isfield(item,'compensate') && item.compensate
                    error('hjb:jump','Compensated Levy jumps not implemented in v2');
                end
                lam=hjb_column(item.rate,N,'jump.rate');
                if isfield(item,'weight'), lam=lam.*hjb_column(item.weight,N,'jump.weight'); end
                if any(lam<0), error('hjb:jump','Jump intensities must be nonnegative'); end
                target=G.S; if isfield(item,'target'), target=item.target; end
                to=reg; if isfield(item,'to'), to=item.to; end
                if ~isscalar(to) || to<1 || to>R || to~=round(to)
                    error('hjb:jump','Invalid destination regime');
                end
                if ~isequal(size(target),[N,d]), error('hjb:jump','Target must be N-by-d'); end
                if strcmp(G.method,'fd')
                    P=hjb_fd_interp(G,target,opts.jump_outside);
                else
                    lower=G.nodes; upper=G.nodes;
                    for ii=1:d
                        if any(target(:,ii)<lower{ii}(1)-1e-12 | target(:,ii)>upper{ii}(end)+1e-12)
                            if strcmp(opts.jump_outside,'clip')
                                target(:,ii)=min(max(target(:,ii),lower{ii}(1)),upper{ii}(end));
                            else
                                error('hjb:jump','Jump outside collocation domain');
                            end
                        end
                    end
                    P=funbas(G.fspace,target,zeros(1,d));
                end
                zcols=(to-1)*N+(1:N);
                AA(:,zcols)=AA(:,zcols)+spdiags(lam,0,N,N)*P;
                AA(:,jcols)=AA(:,jcols)-spdiags(lam,0,N,N)*G.P0;
                if isfield(item,'reward')
                    fi=fi+lam.*hjb_column(item.reward,N,'jump.reward');
                end
            end
        end
        H=fi+AA*state;
        if isfield(model,'adversary') && ~isempty(model.adversary)
            % Nature cannot choose a dynamically inadmissible alternative.
            H(~feasible)=inf;
        else
            H(~feasible)=-inf;
        end
        scores(:,t)=H;
        Llist{t}=AA; Flist{t}=fi; Xlist{t}=xp;
    end
    worst=-inf(N,ncontroller);
    worstIndex=zeros(N,ncontroller);
    for ctrl=1:ncontroller
        ix=find(group==ctrl);
        [worst(:,ctrl),at]=min(scores(:,ix),[],2);
        worstIndex(:,ctrl)=ix(at);
    end
    worst(isinf(worst) & worst>0)=-inf; % no admissible nature response
    [score,bestController]=max(worst,[],2);
    which=worstIndex(sub2ind(size(worstIndex),(1:N)',bestController));
    if any(~isfinite(score))
        ix=find(~isfinite(score),1);
        error('hjb:infeasible','No feasible policy in regime %d at grid node %d',reg,ix);
    end
    p=size(Xlist{1},2);
    if any(cellfun(@(x)size(x,2)~=p,Xlist))
        error('hjb:policy','Different candidates must have same number of controls');
    end
    policy=zeros(N,p); Rop=sparse(N,N*R); fr=zeros(N,1);
    for t=1:nc
        chosen=find(which==t);
        if isempty(chosen), continue; end
        policy(chosen,:)=Xlist{t}(chosen,:);
        Rop(chosen,:)=Llist{t}(chosen,:);
        fr(chosen)=Flist{t}(chosen);
    end
    rows=(reg-1)*N+(1:N);
    A(rows,:)=Rop;
    fval(:,reg)=fr; pol{reg}=policy; branch(:,reg)=which;
    candscores{reg}=worst;
    % Nature's realized policy at each grid node (if provided).
    advwidth=0;
    for t=1:nc
        ax=adversary_x{t};
        if ~isempty(ax)
            if isscalar(ax), ax=repmat(ax,N,1); end
            if size(ax,1)~=N, error('hjb:adversary','Nature x must have N rows'); end
            adversary_x{t}=ax;
            advwidth=max(advwidth,size(ax,2));
        end
    end
    if advwidth>0
        chosenAdv=nan(N,advwidth);
        for t=1:nc
            ix=find(which==t);
            if isempty(ix) || isempty(adversary_x{t}), continue; end
            ax=adversary_x{t};
            if size(ax,2)~=advwidth
                error('hjb:adversary','All nature x must use the same number of columns');
            end
            chosenAdv(ix,:)=ax(ix,:);
        end
        advpol{reg}=chosenAdv;
    else
        advpol{reg}=[];
    end
    advbranch(:,reg)=reshape(adversary_index(which),N,1);
end
out=struct('A',A,'f',fval,'rho',rho,'policy',{pol},'branch',branch,...
    'candidates_count',candidates_count,'max_candidate_scores',{candscores},...
    'adversary',{advpol},'adversary_branch',advbranch);
if strcmp(G.method,'fd') && opts.check_generator
    hjb_check_generator(A,opts.generator_tol,true);
end
end

function col=hjb_column(v,N,name)
if isscalar(v), col=repmat(v,N,1); else, col=v(:); end
if numel(col)~=N || any(~isfinite(col))
    error('hjb:shape','%s must have N finite values',name);
end
end

function q=hjb_policy_q(model,cand,S,x,reg,R,params)
N=size(S,1); q=zeros(N,R);
if isfield(cand,'q') && ~isempty(cand.q)
    v=cand.q;
elseif isfield(model,'Q') && ~isempty(model.Q)
    if isa(model.Q,'function_handle'), v=model.Q(S,x,reg,params{:});
    else, v=model.Q; end
else
    return
end
if ismatrix(v) && isequal(size(v),[R,R])
    q=repmat(v(reg,:),N,1);
elseif ndims(v)==3 && isequal(size(v),[N,R,R])
    for dest=1:R, q(:,dest)=v(:,reg,dest); end
elseif isequal(size(v),[N,R])
    q=v;
elseif isscalar(v) && R==1
    q=zeros(N,1);
else
    error('hjb:Q','Q must be R-by-R, N-by-R-by-R, or N-by-R row');
end
q(:,reg)=0;
if any(~isfinite(q(:))) || any(q(:)<-1e-12)
    error('hjb:Q','Offdiagonal transition rates must be nonnegative');
end
end
