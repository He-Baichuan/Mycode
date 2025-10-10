function [hw, hP, hX, hlambda, hwel] = eqm_iter(ttar, htau, m)
% Initial values
hw = ones(m.N, 1);
hX = ones(m.N, 1);

% iteration parameters
tol = 1e-10;
weight = 0.1;
crit = 0;
diff = 1;

    while diff > tol && crit < 1000
        htar = (1 + ttar)./(1 + m.tar);
        hkappa = htar .* htau;
        psi = (hkappa).^(-m.theta).*repmat((hw).^(-m.theta), [1,m.N]);
        hlambda = psi./repmat(sum(m.lambda.*psi,1), [m.N, 1]);
        hw_out = sum(hlambda.*m.lambda./(1 + ttar).*repmat(hX'.*m.X',[m.N,1]),2)./m.wL;
        hX_out = (hw.*m.wL + sum(ttar./(1+ttar).*hlambda.*m.lambda.*repmat(hX'.*m.X',[m.N, 1]),1)')./m.X;
        
        r = [hw; hX];
        rr = [hw_out; hX_out];
        diff = max(abs(r-rr));
        disp(diff)
        
        wnum = hw(1);%把一个国家的工资支付作为计价品
        hw = hw./wnum;
        hX = hX./wnum;
        
        hw = weight*hw_out + (1 - weight)*hw;
        hX = weight*hX_out + (1 - weight)*hX;
        crit = crit+1;
    end

hP = (sum(m.lambda.*psi,1)').^(-1/m.theta);
hP = hP./wnum;


hwel = (hw.*m.wL + sum(ttar./(1+ttar).*hlambda.*m.lambda.*repmat(hX'.*m.X',[m.N, 1]),1)')./(m.wL + sum(m.tar./(1+m.tar).*m.lambda.*repmat(m.X',[m.N, 1]),1)')./hP;

end