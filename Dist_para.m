mMomentHisgram = zeros(param.nMeasure, param.nS);
gridAMoments = zeros(param.nAQuadrature, param.nMeasure, param.nS);

for j = 1:param.nS
    mMomentHisgram(1,j) = sum(param.gridA.*(dsn(:,j)/sum(dsn(:,j))));
    gridAMoments(:,1,j) = param.gridAQuadrature - mMomentHisgram(1,j);

    for i = 2:param.nMeasure 
        mMomentHisgram(i,j) = sum(((param.gridA - mMomentHisgram(1,j)).^i).* (dsn(:,j)/sum(dsn(:,j)) ));
        gridAMoments(:,i,j) = (param.gridAQuadrature - mMomentHisgram(1,j)).^i -mMomentHisgram(i,j);
    end
end
% mMomentHisgram = zeros(param.nMeasure, param.nS);  % 保存每个状态下的矩
% gridAMoments = zeros(param.nAQuadrature, param.nMeasure, param.nS);  % 存储中心化矩
% 
% for j = 1:param.nS
%     % 计算一阶矩：E[a | s_j]
%     mass_j = sum(dsn(:, j));  % 状态 j 下的总质量
%     prob_j = dsn(:, j) / mass_j;  % 归一化后的分布 P(a | s_j)
% 
%     mMomentHisgram(1,j) = sum(param.gridA .* prob_j);  % 一阶矩
%     gridAMoments(:, 1, j) = param.gridAQuadrature - mMomentHisgram(1, j);  % 中心化
% 
%     % 计算高阶矩
%     for i = 2:param.nMeasure
%         mMomentHisgram(i,j) = sum(((param.gridA - mMomentHisgram(1,j)).^i) .* prob_j);  % 高阶矩
%         gridAMoments(:, i, j) = (param.gridAQuadrature - mMomentHisgram(1,j)).^i - mMomentHisgram(i,j);  % 高阶中心矩
%     end
% end


% mHatHistogram = dsn(1,1:end)./sum(dsn,1);
% 对每个状态计算借贷约束（例如 a_min = -phi）下的质量 mHat
mHatHistogram = dsn(1,:) ./ sum(dsn, 1);  % 每列表示每个状态下在约束点 a_min 处的质量

% Solve for market clearing capital stock
f = @(capital) computeMCResidualPolynomials(capital,param,VarAtTime,mMomentHisgram,gridAMoments,mHatHistogram);
 % options = optimoptions('fsolve','Display','off','TolFun',1e-2);
 options = optimset('TolX', 1e-6, 'Display', 'off'); 
if abs(f(eqm_K)) > 1e-4
	[aggregateCapital,err,exitflag] = fzero(f,eqm_K,options);
end
[~,mCoefficients,mParameters,mMoments,mHat] = ...
    computeMCResidualPolynomials(aggregateCapital,mMomentsHistogram,aGridMoments,mHatHistogram);
