function stat=hjb_check_generator(A,tol,throw_on_fail)
%HJB_CHECK_GENERATOR Check conservativity and offdiagonal positivity.
if nargin<2 || isempty(tol), tol=1e-9; end
if nargin<3, throw_on_fail=true; end
n=size(A,1);
rowerr=max(abs(full(sum(A,2))));
off=A-spdiags(diag(A),0,n,n);
z=nonzeros(off);
if isempty(z), mine=0; else, mine=min(z); end
valid=all(isfinite(nonzeros(A))) && rowerr<=tol && mine>=-tol;
stat=struct('valid',valid,'rowSumError',rowerr,'minOffDiagonal',mine);
if ~valid && throw_on_fail
    error('hjb:generator','Invalid Markov generator: row error %.3g, min offdiag %.3g',rowerr,mine);
end
end
