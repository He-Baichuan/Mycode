function [sigyJ,sigxJ]=mom(gx,hx,varshock,J, method)
%[sigyJ,sigxJ]=mom(gx,hx,varshock,J, method)
% 计算无条件方差-协方差矩阵 of x(t) with x(t+J),即 sigxJ=E[x(t)*x(t+J)'], 
% of y(t) with y(t+J),即 sigyJ=E[y(t)*y(t+J)']
%  x(t) 的状态转移方程为
% x(t+1) = hx x(t) + e(t+1)
% y(t) 的政策函数为
% y(t) = gx x(t)
% 这里 Ee(t)e(t)'=varshock
% J是滞后或超前期数
% method =1 : use doubling algorithm
% method neq 1 : use algebraic method


if nargin<4
J=0;
end


if nargin<5
    method =1;
end


if method == 1 
disp('method=doubling')

%Doubling algorithm
hx_old=hx;
sig_old=varshock;
sigx_old=eye(size(hx));
diferenz=.1;
while diferenz>1e-25
sigx=hx_old*sigx_old*hx_old'+sig_old;

diferenz = max(max(abs(sigx-sigx_old)));
sig_old=hx_old*sig_old*hx_old'+sig_old;
hx_old=hx_old*hx_old;
sigx_old=sigx;
end    %不动点迭代


else


%Algebraic method
%计算 x的方差协方差矩阵
disp('method=kronecker')


sigx = zeros(size(hx));
F = kron(hx,hx);
sigx(:) = (eye(size(F))-F)\varshock(:);


end   %if method



%Get E{x(t)*x(t+J)'}
sigxJ=hx^(-min(0,J))*sigx*(hx')^(max(0,J));
% 注意这里区分了是滞后还是超前


%Get E{y(t)*y(t+J)'}
sigyJ=real(gx*sigxJ*gx');
sigxJ=real(sigxJ);