clear; clc; close all;

%% 定义系统参数
m = 10;             
c = 750;            
k = 3.55e5;         

% 基底激励方程参数
a = 100;            
b = 350;            
c_base = 1.96e4;    
F_base = 350;       
omega_base = 100 * pi; 

%% 定义微分方程
odefun = @(t, state) [
    state(2); % dy/dt = y_dot
    (F_base * sin(omega_base * t) - b * state(2) - c_base * state(1)) / a; 
    state(4); % dx/dt = x_dot
    (c * state(2) + k * state(1) - c * state(4) - k * state(3)) / m; 
];

%% 设置初始条件和时间范围
initial_conditions = [0; 0; 0; 0];
t_span = [0, 1];

%% 求解微分方程
options = odeset('RelTol',1e-8,'AbsTol',1e-10); 
[t, state] = ode45(odefun, t_span, initial_conditions, options);

y = state(:,1);
y_dot = state(:,2);
x = state(:,3);
x_dot = state(:,4);

%% 计算传递给基底的力
f_T = c * (y_dot-x_dot) + k * (y-x);

%% 绘制稳态响应和传递力
figure(1);

subplot(2,1,1);
plot(t, x, 'b', 'LineWidth', 1.5);
xlabel('时间 t (秒)', 'FontSize', 12);
ylabel('稳态响应 x(t) (米)', 'FontSize', 12);
title('单自由度系统的稳态响应 x(t)', 'FontSize', 14);
grid on;
xlim([0 1]);

subplot(2,1,2);
plot(t, f_T, 'r', 'LineWidth', 1.5);
xlabel('时间 t (秒)', 'FontSize', 12);
ylabel('传递力 f_T(t) (牛顿)', 'FontSize', 12);
title('传递给基底的力 f_T(t)', 'FontSize', 14);
grid on;
xlim([0 1]);

%% 定义频率范围
omega = linspace(0, 400, 1000); 

%% 计算传递函数
H_yx = (k + 1j*c*omega) ./ (k - m*omega.^2 + 1j*c*omega);

%% 计算幅值和相位
magnitude = abs(H_yx);           
phase = angle(H_yx) * (180/pi);  

%% 绘制幅频特性和相频特性
figure(2);

subplot(2,1,1);
plot(omega, magnitude, 'b', 'LineWidth', 1.5);
grid on;
xlabel('\omega (rad/s)', 'FontSize', 12);
ylabel('|H_{y,x}(\omega)| (m/m)', 'FontSize', 12);
title('传递函数 H_{y,x}(\omega) 的幅频特性', 'FontSize', 14);
xlim([0 400]);

subplot(2,1,2);
plot(omega, phase, 'r', 'LineWidth', 1.5);
grid on;
xlabel('\omega (rad/s)', 'FontSize', 12);
ylabel('\phi(\omega) (°)', 'FontSize', 12);
title('传递函数 H_{y,x}(\omega) 的相频特性', 'FontSize', 14);
xlim([0 400]);
