function L=hjb_fd_operator(G,mu,sigma,boundary,opts)
%HJB_FD_OPERATOR Conservative monotone Cartesian-grid CTMC approximation.
% Supports nonuniform 1D/multidimensional grids with diagonal covariance;
% nonzero cross covariances currently require UNIFORM grids and a local
% positive-type condition (otherwise gives an explicit error).
% boundary: cell(d,1): 'state_constraint' or 'reflecting' per dimension.
N=G.N; d=G.d;
if isempty(mu), mu=zeros(N,d); end
if ~isequal(size(mu),[N,d]) || any(~isfinite(mu(:)))
    error('hjb:drift','Drift must be finite N-by-d');
end
if nargin<4 || isempty(boundary)
    boundary=repmat({'state_constraint'},1,d);
end
if ischar(boundary), boundary=repmat({boundary},1,d); end
if numel(boundary)~=d, error('hjb:boundary','Specify one boundary type per dimension'); end
if nargin<5, opts=hjb_options([],'fd'); end
if isempty(sigma)
    C=zeros(N,d,d);
else
    if d==1 && isvector(sigma), sigma=reshape(sigma,N,1,1); end
    if ismatrix(sigma) && d>1 && size(sigma,2)==d
        sigma=reshape(sigma,N,d,1);
    end
    if size(sigma,1)~=N || size(sigma,2)~=d || any(~isfinite(sigma(:)))
        error('hjb:diffusion','sigma must be finite N-by-d-by-q');
    end
    q=size(sigma,3);
    C=zeros(N,d,d);
    for k=1:d
        for l=k:d
            c=sum(reshape(sigma(:,k,:),N,q).*reshape(sigma(:,l,:),N,q),2);
            C(:,k,l)=c; C(:,l,k)=c;
        end
    end
end
% Build CTMC rates. Self transitions omitted; diagonal set to minus sum.
rows=cell(1,2*d+2*d*(d-1)); cols=rows; vals=rows;
count=0; diagRates=zeros(N,1);
steps=zeros(1,d); isuniform=true;
for k=1:d
    dx=diff(G.nodes{k}); steps(k)=dx(1);
    isuniform=isuniform && max(abs(dx-dx(1)))<=1e-12*max(1,steps(k));
end
cross=zeros(N,d); % variance already allocated to mixed-direction jumps
for k=1:d
    for l=k+1:d
        c=C(:,k,l);
        if max(abs(c))<=1e-13, continue; end
        if ~isuniform
            error('hjb:cross','Mixed diffusion requires uniform grid; use diagonalized state variables');
        end
        r=abs(c)/(2*steps(k)*steps(l));
        cross(:,k)=cross(:,k)+abs(c)*steps(k)/steps(l);
        cross(:,l)=cross(:,l)+abs(c)*steps(l)/steps(k);
        for sg=[-1 1]
            % Positive covariance uses (++),(--); negative uses (+-),(-+).
            sign_l=sg.*sign(c);
            valid=G.subs(:,k)+sg>=1 & G.subs(:,k)+sg<=G.n(k) & ...
                  G.subs(:,l)+sign_l>=1 & G.subs(:,l)+sign_l<=G.n(l);
            ii=find(valid & r>0);
            if ~isempty(ii)
                jj=ii + sg*G.stride(k) + sign_l(ii)*G.stride(l);
                count=count+1; rows{count}=ii; cols{count}=jj; vals{count}=r(ii);
                diagRates(ii)=diagRates(ii)+r(ii);
            end
        end
    end
end
for k=1:d
    idx=G.subs(:,k); lower=idx==1; upper=idx==G.n(k);
    % State constraints are meaningful only without nondegenerate normal diffusion.
    if strcmpi(boundary{k},'state_constraint')
        if any(mu(lower,k)<-opts.infeasible_tol) || any(mu(upper,k)>opts.infeasible_tol)
            error('hjb:boundary','Outward drift at state-constraint boundary, dim %d',k);
        end
        if any(C(lower,k,k)>opts.infeasible_tol) || any(C(upper,k,k)>opts.infeasible_tol)
            error('hjb:boundary','Nonzero normal diffusion requires reflecting boundary, dim %d',k);
        end
    elseif ~strcmpi(boundary{k},'reflecting')
        error('hjb:boundary','Unknown boundary type: %s',boundary{k});
    end
    variance=reshape(C(:,k,k),N,1)-cross(:,k);
    if any(variance < -opts.generator_tol)
        error('hjb:covariance','Positive-type mixed stencil fails in dimension %d; modify grid or covariance',k);
    end
    variance=max(variance,0);
    dx=diff(G.nodes{k});
    hf=nan(N,1); hb=nan(N,1);
    fidx=idx<G.n(k); bidx=idx>1;
    hf(fidx)=dx(idx(fidx)); hb(bidx)=dx(idx(bidx)-1);
    fp=zeros(N,1); bp=zeros(N,1);
    inside=fidx & bidx;
    fp(inside)=variance(inside)./(hf(inside).*(hf(inside)+hb(inside)));
    bp(inside)=variance(inside)./(hb(inside).*(hf(inside)+hb(inside)));
    % Neumann reflecting boundary via mirrored ghost nodes.
    fp(lower)=variance(lower)/(dx(1)^2);
    bp(upper)=variance(upper)/(dx(end)^2);
    fp(fidx)=fp(fidx)+max(mu(fidx,k),0)./hf(fidx);
    bp(bidx)=bp(bidx)+max(-mu(bidx,k),0)./hb(bidx);
    ii=find(fidx & fp>0);
    if ~isempty(ii)
        count=count+1; rows{count}=ii; cols{count}=ii+G.stride(k); vals{count}=fp(ii);
        diagRates(ii)=diagRates(ii)+fp(ii);
    end
    ii=find(bidx & bp>0);
    if ~isempty(ii)
        count=count+1; rows{count}=ii; cols{count}=ii-G.stride(k); vals{count}=bp(ii);
        diagRates(ii)=diagRates(ii)+bp(ii);
    end
end
rr=(1:N)'; rows{count+1}=rr; cols{count+1}=rr; vals{count+1}=-diagRates;
L=sparse(vertcat(rows{1:count+1}),vertcat(cols{1:count+1}),...
         vertcat(vals{1:count+1}),N,N);
end
