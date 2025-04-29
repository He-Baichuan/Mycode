function [Vyr,Vxr,Vy,Vx]=variance_decomposition(gx,hx,ETA1,ETA2)
%[Vy,Vx]=variance_decomposition(gx,hx,ETA1)
% 计算 y 和 x的方差分解, 
%
% x_t+1 = hx x_t + ETA1 epsilon_t+1
% y_t = gx x_t
% epsilon_t ~ iidN(0,I)
%
%[Vy,Vx] = variance_decomposition(gx,hx,ETA1,ETA2) 
%
% x_t+1 = hx x_t + ETA1 epsilon_t+1
% y_t = gx x_t + ETA2 mu_t
% epsilon_t ~ iidN(0,I) and mu_t ~ iidN(0,I) 
% E(epsilon_t mu_T) = 0 for all t,T


Vy = [];
Vx = [];


n1 = size(ETA1,2);

for j=1:n1
I1 = zeros(n1);

I1(j,j) = 1;

V1 = ETA1*I1*ETA1';

[sigy,sigx]=mom(gx,hx,V1);

Vy = [Vy;diag(sigy)'];
Vx = [Vx;diag(sigx)'];

end

if nargin==4 %否则不包含政策函数的冲击

n2 = size(ETA2,2);

for j=1:n2
I2 = zeros(n2);

I2(j,j) = 1;

V2 = ETA2*I2*ETA2';

Vy = [Vy;diag(V2)'];

end %for j=1:n2

end %if nargin==4

Vyr = Vy*diag(sum(Vy))^(-1);
Vxr = Vx*diag(sum(Vx))^(-1);