clc
clear

load test_tc.mat tc;
load window_size_results.mat window_size;

j = 4; k = 3;
dFC_4_3 = compute_dfc_from_window_size(tc, window_size, j, k, 'MinL', 7, 'Method', 'weighted');
plot(dFC_3_4);
xlabel('Time (TR)');
ylabel('dFC (corr)');
title('Dynamic FC between channel 4 and 3');
