function [h1, h2, h3, h4] = function_simulation_Zhuang_2019_NI_separate(f1, f2, f3, TR, tf_sec, fixed_window_size_seconds)
%% 生成模拟信号并分别绘制四个图

% 设置默认参数（如果输入参数不足）
if nargin < 6
    fixed_window_size_seconds = [30, 60, 120]; % 默认窗口大小
end
if nargin < 5
    tf_sec = 400; % 默认总时长
end
if nargin < 4
    TR = 2; % 默认TR
end
if nargin < 3
    f3 = 1/20; % 默认频率 0.05 Hz
end
if nargin < 2
    f2 = 1/25; % 默认频率 0.04 Hz
end
if nargin < 1
    f1 = 1/50; % 默认频率 0.02 Hz
end

lw = 2;  % lineWidth in plot  
fontSize = 16;
sliding_window_size_vec = ceil(fixed_window_size_seconds/TR);
N_win = numel(sliding_window_size_vec);
t = [1:TR:tf_sec]'; tdim = size(t,1);

%% 生成信号
x1 = zeros(tdim,1);
index1 = 1:floor(200/TR);
x1(index1) = cos(2*pi*f1*t(index1));
index2 = floor(200/TR)+1:tdim;
x1(index2) = cos(2*pi*f3*t(index2));
x2 = cos(2*pi*f2*t);
[x1] = function_energy_normalize(x1);
[x2] = function_energy_normalize(x2);

static_corr = corr(x1,x2);

%% 获取屏幕尺寸以便合理布局
screen_size = get(0, 'ScreenSize');
screen_width = screen_size(3);
screen_height = screen_size(4);
fig_width = 800;
fig_height = 600;

%% 图1: 时间序列 - 左上角
h1 = figure(1); 
set(h1,'position',[50, screen_height-fig_height-50, fig_width, fig_height], 'name', 'Time Series');
plot(t, x1, 'b-', t, x2, 'r-', 'LineWidth', lw);  
grid on;
legend({'y^1 (频率切换信号)', 'y^2 (恒定频率信号)'}, 'fontSize', 12);
title(['时间序列: TR=', num2str(TR), 's, f1=', num2str(f1), 'Hz→f3=', num2str(f3), 'Hz@200s, f2=', num2str(f2), 'Hz'], 'FontSize', fontSize);
xlabel('时间 (秒)', 'FontSize', fontSize);
ylabel('信号强度', 'fontSize', fontSize);
set(gcf, 'color', 'w');

%% 图2: 频率域 - 右上角
h2 = figure(2);
set(h2,'position',[screen_width-fig_width-50, screen_height-fig_height-50, fig_width, fig_height], 'name', 'Frequency Domain');
[Px] = fourier_transform_plot_f_vs_amplitude(x1, TR, tdim);
[Py, f] = fourier_transform_plot_f_vs_amplitude(x2, TR, tdim);
plot(f, Px, 'b-', f, Py, 'r-', 'LineWidth', lw); 
grid on;
legend({'y^1', 'y^2'}, 'fontSize', 12);
title('频域分析', 'FontSize', fontSize);
xlabel('频率 (Hz)', 'FontSize', fontSize);
ylabel('振幅', 'fontSize', fontSize);
set(gcf, 'color', 'w');

%% 计算瞬时频率和周期
x11 = [flip(x1); x1; flip(x1)];
x22 = [flip(x2); x2; flip(x2)];
[x11] = function_energy_normalize(x11);
[x22] = function_energy_normalize(x22);
fx = function_calculate_instantaneous_frequency(x11, 1);
pe1 = 1./(fx(tdim+1:tdim*2));
fy = function_calculate_instantaneous_frequency(x22, 1);
pe2 = 1./fy(tdim+1:tdim*2);
td = round(max([pe1'; pe2'])); td = td(:);

%% 图3: 瞬时周期 - 左下角
h3 = figure(3);
set(h3,'position',[50, 50, fig_width, fig_height], 'name', 'Instantaneous Period');
plot(t, pe1, 'b-', t, pe2, 'r-', 'LineWidth', lw); 
hold on; 
plot(t, td, 'g--', 'LineWidth', lw+1);
grid on;

% 绘制固定窗口大小
legend_display = cell(N_win, 1);
for j = 1:N_win
    plot(t(1:tdim-1), sliding_window_size_vec(j)*ones(tdim-1,1), '--', 'LineWidth', lw);
    legend_display{j,1} = ['固定窗口:', num2str(fixed_window_size_seconds(j)), 's'];
end
hold off;

xlabel('时间 (秒)', 'fontSize', fontSize);
ylabel('周期 (TR)', 'fontsize', fontSize);
legend(['y^1瞬时周期'; 'y^2瞬时周期'; '最优窗口大小'; legend_display], 'fontSize', 10, 'Location', 'best');
title('瞬时周期与窗口大小比较', 'FontSize', fontSize);
set(gcf, 'color', 'w');

%% 图4: 动态功能连接 - 右下角
h4 = figure(4);
set(h4,'position',[screen_width-fig_width-50, 50, fig_width, fig_height], 'name', 'Dynamic Functional Connectivity');

% 计算最优窗口相关性
opt_window_corr = NaN(tdim,1);
for i = 1:tdim-1
    if round(i-td(i,1)/2) >= 1 && round(i+td(i,1)/2) <= tdim
        opt_window_corr(i,1) = corr(x1(round(i-td(i,1)/2):round(i+td(i,1)/2)), x2(round(i-td(i,1)/2):round(i+td(i,1)/2)));
    end
end

% 计算固定窗口相关性
sliding_window_corr = NaN(tdim, N_win);
for j = 1:N_win
    wn = sliding_window_size_vec(j);
    for i = 1:tdim
        if round(i-wn/2) >= 1 && round(i+wn/2) <= tdim
            sliding_window_corr(i,j) = corr(x1(round(i-wn/2):round(i+wn/2)), x2(round(i-wn/2):round(i+wn/2)));
        end
    end
end

% 绘制所有相关性曲线
plot(t, opt_window_corr, 'g-', 'LineWidth', lw); 
hold on;
for j = 1:N_win
    plot(t, sliding_window_corr(:,j), '--', 'LineWidth', lw);
end
plot(t, static_corr*ones(tdim,1), 'm--', 'LineWidth', 1.25);
hold off;
grid on;

legend_labels = ['最优自适应窗口'; legend_display; ['静态相关性: ' num2str(static_corr, '%.2f')]];
legend(legend_labels, 'fontSize', 10, 'Location', 'best');
title('动态功能连接 (相关性)', 'FontSize', fontSize);
xlabel('时间 (秒)', 'fontSize', fontSize);
ylabel('皮尔逊相关系数', 'fontSize', fontSize);
set(gcf, 'color', 'w');

end

%% 辅助函数（保持不变）
function [TCEN,energy] = function_energy_normalize(TC)
TC = TC(:);
tdim = size(TC,1);
TCEN = TC./sqrt(sum(TC.^2)./tdim);
energy = 1/tdim*sum(TCEN.^2);
end

function [f,period] = function_calculate_instantaneous_frequency(x1,TR)
P1_hilbert = hilbert(x1);
tdim = numel(x1);
f = abs(1/TR/(2*pi)*diff(unwrap(angle(P1_hilbert))));
period = 1./f;
period(period>(TR*tdim)) = 60;
end

function [P1,f] = fourier_transform_plot_f_vs_amplitude(ft,deltat,L)
ft = ft(:);
if mean(ft) ~= 0
    ft = ft - mean(ft);
end
Fs = 1/deltat;
t = (0:L-1)'*deltat;
Fw = fft(ft,L);
P2 = abs(Fw/L);
P1 = P2(1:round(L/2));
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*((0:(ceil(L/2)-1))')/L;
end