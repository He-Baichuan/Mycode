function [gx,hx,exitflag]=gx_hx(fy,fx,fyp,fxp,stake)
%[gx,hx,exitflag]=gx_hx(fy,fx,fyp,fxp,stake);

%计算一阶扰动的雅可比矩阵 hx gx
%基于Schmitt-Grohe and Uribe (JEDC, 2004)
%格式：
%E_t[f(yp,y,xp,x)=0.
%状态空间表示的解（基于隐函数定理）
%xp = h(x,sigma) + sigma * eta * ep
%y = g(x,sigma).
%在稳态点 (x,sigma)=(xbar,0)附近， 近似政策函数和状态转移方程 g  h，
%  
%h(x,sigma) = xbar + hx *(x-xbar) 其中 xbar=h(xbar,0),
%
%g(x,sigma) = ybar + gx * (x-xbar),
%其中 ybar=g(xbar,0). 

% exitflag： 取值 0 (no solution), 1 (unique solution), 2 (indeterminacy), or 3 (z11 is not invertible).
% Inputs: fy fyp fx fxp stake
% stake： 设定了 hx 的所有特征值模长的上界(其中默认值 stake=1).
% Outputs: gx hx exitflag

tol = 1e-10;
if nargin<5
    stake=1;
end
exitflag = 1;

%构造系统矩阵 A,B
A = [-fxp -fyp];
B = [fx fy];
NK = size(fx,2);%非前定变量的个数

%广义schur分解
[s,t,q,z] = qz(A,B);   
   % s(abs(s) < tol) = 0;
   % t(abs(t) < tol) = 0;
   % q(abs(q) < tol) = 0;
   % z(abs(z) < tol) = 0;

%取平稳根 
sroot = (abs(diag(t))<stake*abs(diag(s)));  
nk=sum(sroot);%平稳根的个数

%将平稳根重排在左上角
[s,t,q,z] = ordqz(s,t,q,z,sroot);   

%矩阵分块
z21 = z(nk+1:end,1:nk);
z11 = z(1:nk,1:nk);

s11 = s(1:nk,1:nk);
t11 = t(1:nk,1:nk);

%识别BK条件
if nk>NK
    warning('The Equilibrium is Locally Indeterminate')
    exitflag=2;
elseif nk<NK
    warning('No Local Equilibrium Exists')
    exitflag = 0;
end

if rank(z11)<nk;
    warning('Invertibility condition violated')
    exitflag = 3;
end

%计算求解
z11i = z11\eye(nk);
gx = real(z21*z11i);  
hx = real(z11*(s11\t11)*z11i);
