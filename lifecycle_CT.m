clear; close all; clc;

参数设置
tic
ga = 2;
sig2 = 0.8^2;
corr = exp(-0.9);
rho = 0.05;
r = 0.035;
w = 1;

z_mean = exp(sig2/2);
J = 15;
zmax = 2.5;
zmin = 0.75;
amin = 0;
amax = 100;
I = 300;

T = 75;
N = 300;
dt = T/N;
maxit = 100;
crit = 1e-12;
theta = -log(corr);
a = linspace(amin,amax,I)';
da = (amax-amin)/(I-1);
z = linspace(zmin,zmax,J);
dz = (zmax-zmin)/(J-1);
dz2 = dz^2;
aa = repmat(a,[1,J]);
zz = repmat(z,[I,1]);
mu = (-theta.*log(z)+sig2/2).*z;
s2 = sig2.*z.^2;
%前向差分（对于状态变量k）
Vaf = zeros(I,J);             
Vab = zeros(I,J);
c = zeros(I,J);
%构造状态变量z的状态转移方程各项的系数矩阵
yy = - s2/dz2 - mu/dz;
chi =  s2/(2*dz2);
zeta = mu/dz + s2/(2*dz2);

%This will be the upperdiagonal of the matrix Aswitch
updiag=zeros(I,1); %This is necessary because of the peculiar way spdiags is defined.
for j=1:J
    updiag=[updiag;repmat(zeta(j),I,1)];
end

%This will be the center diagonal of the matrix Aswitch
centdiag=repmat(chi(1)+yy(1),I,1);
for j=2:J-1
    centdiag=[centdiag;repmat(yy(j),I,1)];
end
centdiag=[centdiag;repmat(yy(J)+zeta(J),I,1)];

%This will be the lower diagonal of the matrix Aswitch
lowdiag=repmat(chi(2),I,1);
for j=3:J
    lowdiag=[lowdiag;repmat(chi(j),I,1)];
end

%Add up the upper, center, and lower diagonal into a sparse matrix
Aswitch=spdiags(centdiag,0,I*J,I*J)+spdiags(lowdiag,-I,I*J,I*J)+spdiags(updiag,I,I*J,I*J);

%Preallocation
v = zeros(I,J,N);
% terminal condition on value function: value of death \approx 0
small_number1 = 10^(-8); small_number2 = 10^(-8);
v_terminal = small_number1*(small_number2 + aa).^(1-ga)/(1-ga);

主循环
V = v_terminal;
for n=N:-1:1
        disp(['age = ', num2str(n*dt)])
        v(:,:,n)=V;
        % forward difference
        dVf(1:I-1,:) = (V(2:I,:)-V(1:I-1,:))/da;
        dVf(I,:) = (w*z + r.*amax).^(-ga); %will never be used, but impose state constraint a<=amax just in case
        % backward difference
        dVb(2:I,:) = (V(2:I,:)-V(1:I-1,:))/da;
        dVb(1,:) = (w*z + r.*amin).^(-ga); %state constraint boundary condition
        %consumption and savings with forward difference
        cf = dVf.^(-1/ga);
        ssf = w*zz + r.*aa - cf;
        %consumption and savings with backward difference
        cb = dVb.^(-1/ga);
        ssb = w*zz + r.*aa - cb;
        %consumption and derivative of value function at steady state
        c0 = w*zz + r.*aa;
        %upwind method
        If = ssf > 0; %positive drift --> forward difference
        Ib = ssb < 0; %negative drift --> backward difference
        I0 = (1-If-Ib); %at steady state

        c = cf.*If + cb.*Ib + c0.*I0;
        u = c.^(1-ga)/(1-ga);


        %CONSTRUCT MATRIX
        X = - min(ssb,0)/da;
        Y = - max(ssf,0)/da + min(ssb,0)/da;
        Z = max(ssf,0)/da;


        updiag = 0;
        for j=1:J
            updiag=[updiag;Z(1:I-1,j);0];
        end
        
        centdiag = reshape(Y,I*J,1);
        
        lowdiag=X(2:I,1);
        for j=2:J
            lowdiag=[lowdiag;0;X(2:I,j)];
        end
        
        A=Aswitch+spdiags(centdiag,0,I*J,I*J)+spdiags([updiag;0],1,I*J,I*J)+spdiags([lowdiag;0],-1,I*J,I*J);

       if max(abs(sum(A,2)))>10^(-9)
           disp('Improper Transition Matrix')
           break
       end

        %%Note the syntax for the cell array
        A_t{n} = A;
        B = (1/dt + rho)*speye(I*J) - A;
        
        u_stacked = reshape(u,I*J,1);
        V_stacked = reshape(V,I*J,1);
        
        b = u_stacked + V_stacked/dt;
        V_stacked = B\b; %SOLVE SYSTEM OF EQUATIONS
        
        V = reshape(V_stacked,I,J);

        c_t{n} = c;
        ss_t{n} = w*zz + r.*aa - c;
end
toc
plot(a,c_t{1})
plot(a,ss_t{1},a,zeros(I,1),'k--')

plot(a,c_t{N-1})
plot(a,ss_t{N-1},a,zeros(I,1),'k--')
subplot(1,2,1)
set(gca,'FontSize',16)
plot(a,c_t{1/dt}(:,1),a,c_t{1/dt}(:,J),a,c_t{40/dt}(:,1),a,c_t{40/dt}(:,J),a,c_t{70/dt}(:,1),a,c_t{70/dt}(:,J),'Linewidth',2)
legend('Age 1, Lowest Income','Age 1, Highest Income','Age 40, Lowest Income','Age 40, Highest Income','Age 70, Lowest Income','Age 70, Highest Income')
ylim([0 5])
xlabel('Wealth')
ylabel('Consumption, c(a,z,t)')

subplot(1,2,2)
set(gca,'FontSize',16)
plot(a,ss_t{1/dt}(:,1),a,ss_t{1/dt}(:,J),a,ss_t{40/dt}(:,1),a,ss_t{40/dt}(:,J),a,ss_t{70/dt}(:,1),a,ss_t{70/dt}(:,J),a,zeros(I,1),'k--','Linewidth',2)
legend('Age 1, Lowest Income','Age 1, Highest Income','Age 40, Lowest Income','Age 40, Highest Income','Age 70, Lowest Income','Age 70, Highest Income')
ylim([-5 2])
xlabel('Wealth')
ylabel('Saving, s(a,z,t)')
