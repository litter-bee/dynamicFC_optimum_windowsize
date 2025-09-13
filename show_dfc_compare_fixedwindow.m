clc
clear

load test_tc.mat tc;
load window_size_results.mat window_size;

TR = 0.72;
ts1 = tc(:,3);
ts2 = tc(:,4);
N = length(ts1);

% 固定窗口大小（秒 -> 点数）
win_lengths_sec = [30, 60, 90];
win_lengths_pts = round(win_lengths_sec / TR);

% 保存结果
DFC_fixed = zeros(N, length(win_lengths_pts));

% 计算固定窗口 DFC
for w = 1:length(win_lengths_pts)
    win = win_lengths_pts(w);
    tmp_corr = nan(N,1);
    half_w = floor(win/2);
    for t = 1:N
        idx1 = max(1, t-half_w);
        idx2 = min(N, t+half_w);
        if idx2-idx1+1 >= 5  % 至少保证窗口长度大于5
            r = corr(ts1(idx1:idx2), ts2(idx1:idx2));
            tmp_corr(t) = r;
        end
    end
    DFC_fixed(:,w) = tmp_corr;
end

% 动态窗口
dyn_win = round(window_size{4,3}); % 动态窗口长度 (单位: 点数)
DFC_dynamic = nan(N,1);
for t = 1:N
    win = dyn_win(t);
    half_w = floor(win/2);
    idx1 = max(1, t-half_w);
    idx2 = min(N, t+half_w);
    if idx2-idx1+1 >= 5
        r = corr(ts1(idx1:idx2), ts2(idx1:idx2));
        DFC_dynamic(t) = r;
    end
end

% 绘制对比图
time = (0:N-1)*TR; % 时间轴（秒）
figure; hold on;
% plot(time, DFC_fixed(:,1), 'r-', 'LineWidth',1.5); % 30s
plot(time, DFC_fixed(:,2), 'g-', 'LineWidth',1.5); % 60s
% plot(time, DFC_fixed(:,3), 'b-', 'LineWidth',1.5); % 90s
plot(time, DFC_dynamic, 'k--', 'LineWidth',2);     % 动态
legend('60s Fixed','Dynamic');
% legend('30s Fixed','60s Fixed','90s Fixed','Dynamic');
xlabel('Time (s)'); ylabel('DFC (Correlation)');
title('DFC Comparison for Channel 3 & 4');
grid on;
