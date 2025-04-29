function grids = makeHAGrids(a,na,e,ne, Pi_e)
   % 构造一个 HAGrids 结构体
    % 输入：
    %   a    - 资产网格向量（Nx×1 或 1×Nx）
    %   e    - 冲击状态向量（Mx×1 或 1×Mx）
    %   Pi_e - 冲击状态转移矩阵（Mx×Mx）
    %
    % 输出：
    %   grids - 包含字段 .a, .na, .e, .ne, .Pi_e 的 struct

    grids.a    = a(:);            % 保证列向量
    grids.na   = na;
    grids.e    = e(:);
    grids.ne   = ne;
    grids.Pi_e = Pi_e;
end
