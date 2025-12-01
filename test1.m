clc;
crit = 10;
tol = 1e-6;
iteration = 1;
maxit = 500;
option2 = optimoptions(@fminunc,'Algorithm','quasi-newton','Display','notify-detailed',...
	'MaxFunEvals',50000,'TolFun',1e-12,'GradObj','on','MaxIter',1000);
mMoments = mMomentHisgram;
% gridAMoments = gridAMoments;
mHat = mHatHistogram;
% [ap, ~, ~] = SolveHHprobEGM(param, r, w, tau_l, tau_k, transfer);
ap = decision.ap;
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
    % j = 1;
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

 end



% parametersResidual(zeros(param.nMeasure,1),squeeze(gridAMoments(:,:,1)),param);

