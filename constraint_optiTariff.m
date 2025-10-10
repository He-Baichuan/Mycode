function [c,ceq] = constraint_optiTariff(x,i,m)
% variables
hX = x(1: m.N);
hP = x(m.N + 1: 2*m.N);
hw = x(2*m.N + 1: 3*m.N);
ttar = x(3*m.N + 1: 4*m.N);

% tariff matrix
ttar_mat = m.tar;
ttar_mat(:, i) = ttar;
htar = (1 + ttar_mat)./(1 + m.tar);
psi = htar.^(-m.theta).*repmat((hw).^(-m.theta), [1,m.N]);
hlambda = psi./repmat(sum(m.lambda.*psi,1), [m.N, 1]);

hw_out = sum(hlambda.*m.lambda./(1 +ttar_mat).*repmat(hX'.*m.X',[m.N,1]),2)./m.wL;
hX_out = (hw.*m.wL + sum(ttar_mat./(1+ttar_mat).*hlambda.*m.lambda.*repmat(hX'.*m.X',[m.N, 1]),1)')./m.X;
hP_out = (sum(m.lambda.*psi,1)').^(-1/m.theta);

ceq = [hw - hw_out; hX - hX_out; hP - hP_out; hw(1) - 1; ttar(i)]; 
c = [];
end