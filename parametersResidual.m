% function values = parametersResidual(parameters, gridMoments, param)
function [values, derivations] = parametersResidual(parameters, gridMoments, param)

% copy by Winberry
% 计算最优系数来参数化分布
% 输入：parameters 参数的猜测值 gridMoments 矩条件 param 参数
% 输出：values 函数值 derivations 导数值
% Weights = param.QuadratureWeights;
% nMeasure = param.nMeasure;
% %函数值计算
% values = Weights' * gridMoments *parameters;
% %导数值计算
%     if nargout>1
%         derivations = sum(repmat(Weights,[1 nMeasure]) .* gridMoments .* repmat(exp(gridMoments * parameters),...
% 		[1 nMeasure]),1)';
% 
%     end
    Weights = param.QuadratureWeights;          % nQuad × 1
    z       = gridMoments * parameters;         % nQuad × 1
    ez      = exp(z);                           % nQuad × 1

    % 目标：近似归一化常数 ∫ exp(φ(a)'θ) da
    values = Weights' * ez;                     % 标量

    if nargout > 1
        % 梯度：∂/∂θ  ∫ exp(φ(a)'θ) da
        % ≈ gridMoments' * (Weights .* exp(gridMoments * θ))
        derivations = sum(repmat(Weights,[1 param.nMeasure]) .* gridMoments .* repmat(exp(gridMoments * parameters),...
 		[1 param.nMeasure]),1)';   % nMeasure × 1
    end
end
