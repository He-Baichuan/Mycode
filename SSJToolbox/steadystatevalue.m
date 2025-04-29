function val = steadystatevalue(s, m, endogexog, which)
% STEADYSTATEVALUE  Get steady-state of variable s
%  val = steadystatevalue(s, m, endogexog, which)
    if nargin<4, which='terminal'; end
    if strcmp(endogexog,'endog')
        vec   = m.ss.(which);
        names = m.vars('0');
    elseif strcmp(endogexog,'exog')
        vec   = m.ss.exog;
        names = m.vars('exog');
    else
        error('endogexog must be ''endog'' or ''exog''.');
    end
    idx = find(strcmp(names,s));
    assert(~isempty(idx), 'Variable not found');
    val = vec(idx);
end