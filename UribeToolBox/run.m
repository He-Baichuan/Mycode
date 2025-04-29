clear;
close all;
clc;
%% 参数设定
beta = 0.99;
theta = 0.36;
delta = 0.025;
gamma = 0.95;
A = 2;
h0 = 0.583;
B = A*log(1-h0)/h0;
approx = 1;

%% 模型设定
syms c cp k kp r rp lambda lambdap h hp y yp

f1 = cp-beta*(c*(rp+1-delta));
f2 = c+(1-theta)*y/(B*h);
f3 = c-y-(1-delta)*k+kp;
f4 = y-lambda*k^(theta)*h^(1-theta);
f5 = r-theta*lambda*k^(theta-1)*h^(1-theta);
f6 = log(lambdap)-gamma*log(lambda);

f = [f1;f2;f3;f4;f5;f6];
x = [k,lambda];
y = [c,r,h,y];
yp = [cp,rp,hp,yp];
xp = [kp,lambdap];

f = subs(f, [x,y,xp,yp], (exp([x,y,xp,yp])));%相对偏离程度
[fx,fxp,fy,fyp,fypyp,fypy,fypxp,fypx,fyyp,fyy,fyxp,fyx,fxpyp,fxpy,fxpxp,fxpx,fxyp,fxy,fxxp,fxx] = anal_deriv(f,x,y,xp,yp,approx);

%% 稳态设定
kbar = 12.6698;
hbar = 0.3335;
ybar = 1.2353;
cbar = 0.9186;
rbar = 0.0351;
wbar = 2.3706;
c = log(cbar);
cp = c;
k = log(kbar);
kp = k;
r = log(rbar);
rp = r;
lambda = 0;
lambdap = lambda;
h = log(hbar);
hp = h;
y = log(ybar);
yp = y;
num_deriv;
[gx,hx] = gx_hx(nfy,nfx,nfyp,nfxp);
s0 = [0,1];
T = 100;
[IR,IRy,IRx,tt] = Irf(gx,hx,s0,T);

figure
plot(tt,IRx(1,:),tt,IRx(2,:));
figure
plot(tt,IRy(1,:),tt,IRy(2,:),tt,IRy(3,:),tt,IRy(4,:));
varshock = eye(2,2);
[sigyJ,sigxJ] = mom(gx,hx,varshock,0, 2);

[Y,X, e] = simu_1st(gx, hx,varshock, 800);
figure
plot(Y);
Y = Y(300:end,:); 
sigY = std(Y);
[Vy,Vx]=variance_decomposition(gx,hx,varshock);