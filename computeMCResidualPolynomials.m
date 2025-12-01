function [residual,mCoefficientsOptional,mParametersOptional,mMomentsOptional,mHatOptional] = ...
	computeMCResidualPolynomials(Kguess,param, VarAtTime,mMoments,gridAMoments,mHat)% ExcessDemand
% Computes residual of market-clearing condition, parametric family to approximate distribution
% 
% Inputs
%   (1) capital: candidate aggregate capital stock
%	(2) mMoments: intial guess of moments of distribution (nEpsilon x nMeasure)
%	(3) aGridMoments: grid of centralized moments for computing PDF, corresponding to mMoments
%		(nEpsilon x nAssetsQuadrature x nMoments)
%	(4) mHat: initial guess of mass at borrowing constraint
%
% Outputs
%   (1) residual: residual of market clearing condition
%   (2) (optional) mCoefficientsOptional: coefficients on individual decisions (splines or polynomials)
%   (3) (optional) mParametersOptional: parameters of density away from borrowing constraint
%   (4) (optional) mMomentsOptional: moments of density away from borrowing constraint
%	(5) (optional) mHatOptional: mass at borrowing constraint
% 
    %% 解包参数
    Lbar  = param.Lbar;
    alpha = param.alpha;
    delta = param.delta;
    nA    = param.nA;
    nS    = param.nS;
    gridA = param.gridA(:);      % nA×1

    %% 解包当期宏观变量
    Zbar = VarAtTime(1);   % 总生产率
    tau_l = VarAtTime(2);  % 劳动税率
    tau_k = VarAtTime(3);  % 资本税率
    Gexp  = VarAtTime(4);  % 政府支出占 Y 比例

    %% 1. 给定 Kguess，计算价格与产出
    % r = Zbar * alpha * K^(alpha-1) * L^(1-alpha) - delta
    r = Zbar * alpha * Kguess^(alpha - 1) * Lbar^(1 - alpha) - delta;

    % w = Zbar * (1-alpha) * K^alpha * L^(-alpha)
    w = Zbar * (1 - alpha) * Kguess^alpha * Lbar^(-alpha);

    % Y = Zbar * K^alpha * L^(1-alpha)
    Y = Zbar * Kguess^alpha * Lbar^(1 - alpha);

    % 政府支出 G = Gexp * Y
    G = Gexp * Y;

    %% 2. 政府预算平衡，得到一次性转移 transfer
    % transfer = tau_l * w * Lbar + tau_k * r * Kguess - G
    transfer = tau_l * w * Lbar + tau_k * r * Kguess - G;

    %% 3. 家庭问题：EGM 解出最优决策
    % decision.ap: a'(a,s)，decision.cp: c(a,s) 等
    [ap, ~, ~] = SolveHHprobEGM(param, r, w, tau_l, tau_k, transfer);
    	% Compute weights
	[vIndicesBelow,vIndicesAbove,vWeightBelow,vWeightAbove] = computeLinearWeights(param.gridA,param.gridAQuadrature);
		
	% Linear interpolation
	mAssetsPrimeQuadrature = ap(vIndicesBelow,:) .* repmat(vWeightBelow,1,param.nS) + ...
		ap(vIndicesAbove,:) .* repmat(vWeightAbove,1,param.nS);
    % Compute weights
	[vIndicesBelow,vIndicesAbove,vWeightBelow,vWeightAbove] = computeLinearWeights(param.gridA,-param.phi);
		
	% Linear interpolation
	mAssetsPrimeBC = ap(vIndicesBelow,:) .* repmat(vWeightBelow,1,param.nS) + ...
		ap(vIndicesAbove,:) .* repmat(vWeightAbove,1,param.nS);
    crit = 100; iteration = 1; 
option2 = optimoptions(@fminunc,'Algorithm','quasi-newton','Display','notify-detailed',...
	'MaxFunEvals',50000,'TolFun',1e-12,'GradObj','on','MaxIter',1000);
tol = 1e-8;
maxit = 500;

while crit>tol && iteration<maxit
    mParameters = zeros(param.nMeasure+1,param.nS);
    for j = 1:param.nS
        objParametersResidual = @(x) parametersResidual(x,squeeze(gridAMoments(:,:,j)),param);
        [Parameters, nomralization] = fminunc(objParametersResidual,zeros(param.nMeasure,1),option2);
        mParameters(:,j) = [1/nomralization; Parameters];
    end

%计算 law of motion
    mMomentsNew = zeros(param.nMeasure, param.nS);
    gridAMomentsNew = zeros(param.nAQuadrature, param.nMeasure, param.nS);
%计算一阶矩
    for j = 1:param.nS
        mMomentsNew(1,j) = 0;
        for jp = 1:param.nS

            mMomentsNew(1,j) = mMomentsNew(1,j)+ (1-mHat(1,jp))* param.inv_S(jp)* param.transS(jp,j)*...
                mParameters(1,jp)*param.QuadratureWeights' * (mAssetsPrimeQuadrature(:,jp).*...
                exp(squeeze(gridAMoments(:,:,jp))* mParameters(2:param.nMeasure+1,jp)))+ mHat(1,jp)*...
                param.inv_S(jp)* param.transS(jp,j) * mAssetsPrimeBC(1,jp);
        end
        mMomentsNew(1, j) = mMomentsNew(1,j)/param.inv_S(j);
        gridAMomentsNew(:,1,j) = param.gridAQuadrature - mMomentsNew(1,j);
%计算高阶矩
        for i = 2:param.nMeasure
            mMomentsNew(i,j) = 0;

            for jp = 1:param.nS

                mMomentsNew(i,j) = mMomentsNew(i,j) + (1 - mHat(1,jp)) * param.inv_S(jp) * param.transS(jp,j)*...
                    mParameters(1,jp)*param.QuadratureWeights'*(((mAssetsPrimeQuadrature(:,jp)-...
                    mMomentsNew(1,j)).^i).* exp(squeeze(gridAMoments(:,:,jp))*...
                    mParameters(2:param.nMeasure+1,jp)))+mHat(1,jp)*param.transS(jp,j)*param.inv_S(jp)*...
                    ((mAssetsPrimeBC(1,jp)- mMomentsNew(1,j)).^i);

            end
            mMomentsNew(i,j)= mMomentsNew(i,j)/param.inv_S(j);
            gridAMomentsNew(:,i,j) = (param.gridAQuadrature - mMomentsNew(1,j)).^i - mMomentsNew(i,j);

        end

    end

    %%借贷约束束紧的mass
 mHatNew = zeros(1,param.nS);

    for j = 1:param.nS
        for jp = 1:param.nS
            mHatNew(1,j) = mHatNew(1,j)+ (1 - mHat(1,jp))* param.transS(jp,j)* param.inv_S(jp)*...
                mParameters(1,jp)*param.QuadratureWeights'*((mAssetsPrimeQuadrature(:,jp)<= -param.phi+1e-8).*...
                exp(squeeze(gridAMoments(:,:,jp)) * mParameters(2:param.nMeasure+1 , jp)) )+ ...
                mHat(1,jp)*param.transS(jp,j)*param.inv_S(jp)*(mAssetsPrimeBC(1,jp)<= -param.phi + 1e-8);
        end
        mHatNew(1,j) = mHatNew(1,j)/param.inv_S(j);
    end
%%%%%%%%
%迭代更新
%%%%%%%%

crit  = max([max(abs(mMomentsNew(:) - mMoments(:))), max(abs(mHatNew(:) - mHat(:)))]);
iteration = iteration + 1;
mMoments = mMomentsNew;
mHat = mHatNew;
gridAMoments = gridAMomentsNew;
disp(iteration)
disp(crit)
end


%----------------------------------------------------------------
% Return market clearing residual
%----------------------------------------------------------------

capitalNew = mMoments(1,:) * (param.inv_S .* (1 - mHat)) - param.phi  * ones(1,param.nS)*(param.inv_S .* mHat);
residual = Kguess - capitalNew;

% Also return optional outputs if requested
if nargout > 2

    mParametersOptional = mParameters;
    mMomentsOptional = mMoments;
	mHatOptional = mHat;



end
