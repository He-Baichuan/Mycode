function [IR,IRy,IRx,tt]=Irf(gx,hx,s0,T)
%脉冲响应函数
if nargin < 4
    T =40;
end 
tt = (1:T);
pd = length(s0);
MX=[gx;eye(pd)];
IR=[];
x=s0';
for t=1:T
    IR(:,t)=(MX*x);
    x = hx * x;
end

if nargout>1
    IRx = IR(end-pd+1:end,:);
    IRy = IR(1:end-pd,:);
end


