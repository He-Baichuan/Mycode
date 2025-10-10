function f = obj_optimalTariff(x, i, m)
hX = x(1:m.N);
hP = x(m.N+1:2*m.N);
hw = x(2*m.N + 1: 3*m.N);
ttar = x(3*m.N+1:end);

ttar_mat = m.tar;
ttar_mat(:,i) =ttar;
htar = (1 + ttar_mat)./(1 + m.tar);

psi = htar.^(-m.theta).*repmat((hw).^(-m.theta), [1,m.N]);
hlambda = psi./repmat(sum(m.lambda.*psi,1), [m.N, 1]);
hwel = (hw.*m.wL + sum(ttar_mat./(1 + ttar_mat).*hlambda.*m.lambda.*repmat(hX'.*m.X',[m.N, 1]) ,1 )' )./...
    (m.wL + sum(m.tar./(1 + m.tar).*m.lambda.*repmat(m.X',[m.N, 1]) ,1 )' )./hP; % welfare
f =  -hwel(i);
end
