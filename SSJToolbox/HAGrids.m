function grids = HAGrids( )
% HAGrids  生成资产—冲击网格及其转移矩阵

% 输出 fields:
%   grids.a    - 资产网格向量 (na×1)
%   grids.na   - 网格点个数
%   grids.e    - 劳动 endowment 向量 (ne×1)
%   grids.ne   - 冲击状态个数
%   grids.Pi_e - 转移概率矩阵 (ne×ne)

  

    %——1. 参数设定——
    amin  = 0;        % 借贷约束下限
    Ne    = 2;        % 冲击状态数目
    lamw  = 0.6;      % 找工作概率
    sigma = 0.2;      % 失业概率

    % 转移矩阵：行=当前状态，列=未来状态
    exogTrans = [1-lamw, lamw;
                 sigma,  1-sigma];

    %——2. 劳动 endowment 及归一化——
    e0 = [1.0; 2.5];
    pi_stat = stationarydistribution(exogTrans', 0, 1e5);
    Lbar    = sum(e0 .* pi_stat);
    endow  = e0 ./ Lbar;          % 归一化后加权总劳动=1

    %——3. 资产网格生成——
    agridmin  = amin;
    agridmax  = 400;
    agridsize = 201;
    agrid = linspace(agridmin^(0.25), agridmax^(0.25), agridsize) .^ 4;%凸函数格点
    % agrid = linspace(agridmin, agridmax, agridsize) ;
    %——4. 组装输出——
    grids = makeHAGrids(agrid, agridsize,endow,Ne, exogTrans);
end