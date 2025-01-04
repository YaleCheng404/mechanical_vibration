clear; clc; close all;

% 定义质量矩阵 M
M = [100, 0;
     0, 10];

% 定义刚度矩阵 K
K = [1.96e4 + 3.55e5, -3.55e5;
     -3.55e5, 3.55e5];

%% 求解广义特征值问题 K*v = omega^2*M*v
[Vec, D] = eig(K, M);

% 提取特征值和固有频率
omega_squared = diag(D);
omega = sqrt(omega_squared); % 单位 rad/s

% 转换为频率（Hz）
f = omega / (2*pi);

% 显示固有频率
disp('系统的固有频率（Hz）:');
disp(f);

% 提取振型矩阵
V = Vec;

disp('系统的振型矩阵:');
disp(V);

% 绘制振型图
figure;
subplot(2,1,1);
stem([1, 2], V(:,1), 'filled');
title(['振型 1，频率 = ', num2str(f(1), '%.2f'), ' Hz']);
xlabel('自由度');
ylabel('振型幅值');
xticks([1,2]);
xticklabels({'x_1', 'x_2'});
grid on;

subplot(2,1,2);
stem([1, 2], V(:,2), 'filled');
title(['振型 2，频率 = ', num2str(f(2), '%.2f'), ' Hz']);
xlabel('自由度');
ylabel('振型幅值');
xticks([1,2]);
xticklabels({'x_1', 'x_2'});
grid on;


%% 正则化振型（质量规范化）：v_i^T * M * v_i = 1
V2=V;
for i = 1:size(V2,2)
    norm_factor = sqrt(V2(:,i)' * M * V2(:,i));
    V2(:,i) = V2(:,i) / norm_factor;
end

% 显示正则振型矩阵
disp('正则振型矩阵 V2:');
disp(V2);

% 绘制正则振型图
figure;
for i = 1:size(V2,2)
    subplot(size(V2,2),1,i);
    stem([1, 2], V2(:,i), 'filled');
    title(['正则振型 ', num2str(i), '，频率 = ', num2str(f(i), '%.2f'), ' Hz']);
    xlabel('自由度');
    ylabel(['振型幅值 v_', num2str(i)]);
    xticks([1,2]);
    xticklabels({'x_1', 'x_2'});
    grid on;
end

%% 频率范围设置
f_min = 0;      
f_max = 50;     
num_points = 1000; 

f = linspace(f_min, f_max, num_points); 
omega = 2 * pi * f; 

% 初始化频响函数矩阵
H11 = zeros(1, num_points);
H21 = zeros(1, num_points);
H12 = zeros(1, num_points);

% 计算频响函数矩阵的每个元素
for i = 1:num_points
    w = omega(i);
    Z = K - (w^2) * M; 
    H = inv(Z); 
    H11(i) = H(1,1);
    H21(i) = H(2,1);
    H12(i) = H(1,2);
end

disp('系统的频响函数矩阵:');
disp(H);

% 计算幅值和相位（以度为单位）
H11_mag = abs(H11);
H11_phase = angle(H11) * (180/pi);

H21_mag = abs(H21);
H21_phase = angle(H21) * (180/pi);

H12_mag = abs(H12);
H12_phase = angle(H12) * (180/pi);

% 绘制幅频特性
figure;
subplot(3,2,1);
plot(f, H11_mag, 'b', 'LineWidth', 1.5);
title('|H_{11}(j\omega)|');
xlabel('频率 (Hz)');
ylabel('幅值 |H_{11}|');
grid on;

subplot(3,2,3);
plot(f, H21_mag, 'r', 'LineWidth', 1.5);
title('|H_{21}(j\omega)|');
xlabel('频率 (Hz)');
ylabel('幅值 |H_{21}|');
grid on;

subplot(3,2,5);
plot(f, H12_mag, 'g', 'LineWidth', 1.5);
title('|H_{12}(j\omega)|');
xlabel('频率 (Hz)');
ylabel('幅值 |H_{12}|');
grid on;

% 绘制相频特性
subplot(3,2,2);
plot(f, H11_phase, 'b', 'LineWidth', 1.5);
title('∠H_{11}(j\omega) [度]');
xlabel('频率 (Hz)');
ylabel('相位 ∠H_{11}');
grid on;

subplot(3,2,4);
plot(f, H21_phase, 'r', 'LineWidth', 1.5);
title('∠H_{21}(j\omega) [度]');
xlabel('频率 (Hz)');
ylabel('相位 ∠H_{21}');
grid on;

subplot(3,2,6);
plot(f, H12_phase, 'g', 'LineWidth', 1.5);
title('∠H_{12}(j\omega) [度]');
xlabel('频率 (Hz)');
ylabel('相位 ∠H_{12}');
grid on;

%% 定义固有频率（Hz）和对应的圆频率（rad/s）
f1 = 2.1240;         
f2 = 31.4579;        
omega1 = 2*pi*f1;   
omega2 = 2*pi*f2;   

% 定义激励频率
f_exc = 50;                 
omega = 2*pi*f_exc;        

% 定义外力系数（从正则坐标运动方程得到）
F_modal = [-33.355; -10.605];  

% 计算稳态响应 Q1 和 Q2
Q1 = F_modal(1) / (omega1^2 - omega^2);
Q2 = F_modal(2) / (omega2^2 - omega^2);

% 定义时间向量
t = linspace(0, 0.1, 1000);  

% 计算 q1(t) 和 q2(t)
q1 = Q1 * sin(omega * t);
q2 = Q2 * sin(omega * t);
V2 = [-0.0953, -0.0303;
     -0.0958, 0.3014];

% 计算 x1(t) 和 x2(t)
x1 = V2(1,1)*q1 + V2(1,2)*q2;
x2 = V2(2,1)*q1 + V2(2,2)*q2;

% 计算输入力 f1(t)
f1_t = 350 * sin(omega * t);

% 绘制时间历程图
figure;

% 绘制输入力 f1(t)
subplot(3,1,1);
plot(t, f1_t, 'b', 'LineWidth', 1.5);
title('输入力 f_1(t)');
xlabel('时间 (s)');
ylabel('力幅值 (N)');
grid on;
xlim([0, 0.1]);

% 绘制 x1(t)
subplot(3,1,2);
plot(t, x1, 'r', 'LineWidth', 1.5);
title('稳态响应 x_1(t)');
xlabel('时间 (s)');
ylabel('位移幅值 (m)');
grid on;
xlim([0, 0.1]);

% 绘制 x2(t)
subplot(3,1,3);
plot(t, x2, 'g', 'LineWidth', 1.5);
title('稳态响应 x_2(t)');
xlabel('时间 (s)');
ylabel('位移幅值 (m)');
grid on;
xlim([0, 0.1]);
