clc
clear

load test_tc.mat tc;
load window_size_results.mat window_size;

TR = 0.72;
[N, num_chan] = size(tc);  % N=1186, num_chan=120
ts_ref = tc(:,1);

DFC_matrix = nan(N, num_chan-1);

for j = 2:num_chan
    ts2 = tc(:,j);
    % window_size 是下三角
    if j < 1
        dyn_win = round(window_size{1,j}); 
    else
        dyn_win = round(window_size{j,1});
    end
    
    tmp_corr = nan(N,1);
    for t = 1:N
        win = dyn_win(t);
        half_w = floor(win/2);
        idx1 = max(1, t-half_w);
        idx2 = min(N, t+half_w);
        if idx2-idx1+1 >= 5
            r = corr(ts_ref(idx1:idx2), ts2(idx1:idx2));
            tmp_corr(t) = r;
        end
    end
    DFC_matrix(:, j-1) = tmp_corr;
end

% 绘制热力图
time = (0:N-1)*TR; 
figure;
imagesc(time, 2:num_chan, DFC_matrix');  % 行=通道, 列=时间
set(gca, 'YDir','normal'); 
xlabel('Time (s)');
ylabel('Channel');
title('Dynamic Functional Connectivity (Ref=Channel 1)');
colorbar;
