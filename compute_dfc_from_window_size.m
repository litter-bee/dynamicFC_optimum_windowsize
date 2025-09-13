function dFC_ts = compute_dfc_from_window_size(tc, window_size_cell, j, k, varargin)
% 计算通道 j 与 k 的动态功能连接（dFC）
% tc: tdim x Nq 原始时间序列（已标准化或未标准化均可）
% window_size_cell: Nq x Nq cell 矩阵，元素为长度 tdim 的向量（单位 TR）
% j, k: 通道索引（示例 3,4）
% 可选参数 (Name,Value):
%   'MinL' (default 5) - 最小窗口长度（TR）
%   'Method' (default 'weighted') - 'simple' 或 'weighted'
%   'MakeOdd' (default true) - 窗口长度是否强制为奇数
%
% 返回:
%   dFC_ts: length tdim 的相关系数时间序列

% 解析可选参数
p = inputParser;
addParameter(p,'MinL',5,@(x)isnumeric(x) && x>=1);
addParameter(p,'Method','weighted',@(s)ischar(s));
addParameter(p,'MakeOdd',true,@islogical);
parse(p,varargin{:});
minL = p.Results.MinL;
method = p.Results.Method;
makeOdd = p.Results.MakeOdd;

[tdim, Nq] = size(tc);
if isempty(window_size_cell{j,k})
    error('window_size{%d,%d} is empty.', j, k);
end

ws = window_size_cell{j,k};  % 以 TR 为单位，长度应为 tdim
if numel(ws) ~= tdim
    error('window_size长度与时间点数不一致。');
end

% 把窗口长度从浮点数转换为整数 TR 数
L = round(ws);  % 也可用 ceil/ floor，根据需求
L(L < minL) = minL;  % 强制最小窗长
if makeOdd
    evenIdx = mod(L,2)==0;
    L(evenIdx) = L(evenIdx) + 1;
end

% 取出两通道数据
x_all = tc(:, j);
y_all = tc(:, k);

dFC_ts = nan(tdim,1);

% 主循环：对每个时间点根据 L(t) 计算相关
for t = 1:tdim
    len = L(t);
    half = floor(len/2);
    i1 = t - half;
    i2 = t + half;
    % 边界裁切
    if i1 < 1
        i1 = 1;
    end
    if i2 > tdim
        i2 = tdim;
    end
    x = x_all(i1:i2);
    y = y_all(i1:i2);
    if numel(x) < 3
        dFC_ts(t) = NaN;
        continue;
    end

    if strcmpi(method,'simple')
        % 直接 Pearson（无权重）
        R = corrcoef(x,y);
        r = R(1,2);
        dFC_ts(t) = r;
    else
        % weighted - 使用高斯权重，中心处权重最大
        n = numel(x);
        center = (n+1)/2;
        % sigma 控制窗口权重宽度（经验值）
        sigma = n/6; % 约 99% 在 +/-3sigma
        idx = (1:n);
        w = exp(-0.5*((idx-center).^2)/(sigma^2));
        w = w(:);
        W = sum(w);
        % 加权均值
        mx = sum(w .* x) / W;
        my = sum(w .* y) / W;
        % 加权协方差与方差
        cov_xy = sum(w .* (x - mx) .* (y - my)) / W;
        var_x = sum(w .* (x - mx).^2) / W;
        var_y = sum(w .* (y - my).^2) / W;
        if var_x <= 0 || var_y <= 0
            dFC_ts(t) = NaN;
        else
            dFC_ts(t) = cov_xy / sqrt(var_x * var_y);
        end
    end
end
end
