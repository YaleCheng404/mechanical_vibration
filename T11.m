clear; close all; clc;

%% 参数定义与基本计算
m = 100;                % kg
k = 1.96e4;             % N/m
c = 350;                % N·s/m
F0 = 350;               % 力幅值 N
omega_exc = 100*pi;      % 激励角频率 rad/s

%固有频率与阻尼比
wn = sqrt(k/m);
zeta = c/(2*sqrt(m*k));

fprintf('系统固有频率 wn = %.4f rad/s\n', wn);
fprintf('系统阻尼比 zeta = %.4f\n', zeta);

%% 求稳态响应振幅与相位
%幅值：
H_mag = 1 / sqrt((k - m*omega_exc^2)^2 + (c*omega_exc)^2);
%稳态位移幅值：
Y_amp = F0 * H_mag;

%相位
phi = atan2(-c*omega_exc, (k - m*omega_exc^2));

fprintf('稳态响应幅值约为：Y_amp = %.6e m\n', Y_amp);
fprintf('稳态响应相位(相对激励) phi = %.4f rad\n', phi);

%% 数值求解0~2s内响应
tspan = [0 2];   
y0 = [0;0];      

f = @(t) F0*sin(omega_exc*t);
odefun = @(t,Y)[Y(2); (f(t)-c*Y(2)-k*Y(1))/m];

[t_sol,Y_sol] = ode45(odefun, tspan, y0);
y_sol = Y_sol(:,1);
y_dot_sol = Y_sol(:,2);
y_ddot_sol = diff(y_dot_sol)./diff(t_sol);  
t_accel = t_sol(1:end-1);

f_sol = arrayfun(f, t_sol);

%% 绘图：位移、速度、加速度、激励力时间历程
figure('Name','Time Responses (Force Excitation)','NumberTitle','off');
subplot(4,1,1);
plot(t_sol, f_sol,'LineWidth',1.5);
grid on; xlabel('Time (s)'); ylabel('Force (N)');
title('激励力 f(t)');

subplot(4,1,2);
plot(t_sol, y_sol,'LineWidth',1.5);
grid on; xlabel('Time (s)'); ylabel('Displacement (m)');
title('位移响应 y(t)');

subplot(4,1,3);
plot(t_sol, y_dot_sol,'LineWidth',1.5);
grid on; xlabel('Time (s)'); ylabel('Velocity (m/s)');
title('速度响应 dy/dt(t)');

subplot(4,1,4);
plot(t_accel, y_ddot_sol,'LineWidth',1.5);
grid on; xlabel('Time (s)'); ylabel('Acceleration (m/s^2)');
title('加速度响应 d^2y/dt^2(t)');

%% 频率响应函数幅频与相频特性
omega_vec = linspace(0,300,1000); % rad/s
H_mag_vec = zeros(size(omega_vec));
H_phase_vec = zeros(size(omega_vec));

for i = 1:length(omega_vec)
    w = omega_vec(i);
    H = 1./(k - m*w^2 + 1j*c*w);
    H_mag_vec(i) = abs(H);
    H_phase_vec(i) = angle(H);
end

figure('Name','Frequency Response (Force Excitation)','NumberTitle','off');
subplot(2,1,1);
plot(omega_vec,H_mag_vec,'LineWidth',1.5);
grid on; xlabel('\omega (rad/s)'); ylabel('|H_{y,f}(\omega)| (m/N)');
title('幅频特性');

subplot(2,1,2);
plot(omega_vec,H_phase_vec,'LineWidth',1.5);
grid on; xlabel('\omega (rad/s)'); ylabel('Phase (rad)');
title('相频特性');
