%脉冲响应函数
T = 100;
tt = (1:T);
Rs = size(R,1);
Ps= size(P,1);
shock = zeros(1,T);
shock(1) = 0.01;
lambda_irf = zeros(1,T);
k_irf = zeros(Ps,T);
y_irf = zeros(Rs,T);

for i=1:T
    if i==1
        lambda_irf(:,i) = shock(:,i);
        k_irf(:,i) = Q*lambda_irf(:,i);
        y_irf(:,i) = S*lambda_irf(:,i);
    else
    lambda_irf(:,i) = gamma*lambda_irf(:,i-1)+shock(:,i);
    k_irf(:,i) = P*k_irf(:,i-1) + Q*lambda_irf(:,i);
    y_irf(:,i) = R*k_irf(:,i-1)+S*lambda_irf(:,i);
    end
end



figure
plot(tt,lambda_irf,'-k','LineWidth',1);
plot(tt,zeros(1,T),tt,k_irf,tt,y_irf(1,:),tt,y_irf(2,:),tt,y_irf(3,:),tt,y_irf(4,:),'LineWidth',1)