function [sigy,sigx] = fev(gx,hx,varshock,J)
%[sigy,sigx] = fev(gx,hx,varshock,J)
%计算x(t)、y(t) 在J期后的预测方差矩阵
%sigx=var[x(t+J)-E_tx(t+J)]
%sigy=var(y(t+J)-E_ty(t+J)]
% x(t+1) = hx x(t) + e(t+1)
% y(t) = gx x(t)
% Ee(t)e(t)'=varshock
%J是超前期数



if nargin<3
J=8;
end


sigx = hx*0;

for j=1:J
sigx = hx*sigx*hx' + varshock;
end

sigy = gx*sigx*gx';