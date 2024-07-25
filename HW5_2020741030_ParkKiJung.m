clc, clear
close all

%% Design and Filter
delta_t = 0.01;
sim_time    = 0:delta_t:5;

Pos_pred = 0;
Pos_pred2 = 0;
Vel_pred = 80; 
v_noise = 10 * randn(size(sim_time));     
Vel_measure = 80 + v_noise;

Pos_ori = array_integral(Vel_pred*ones(1, length(sim_time)), delta_t);      %Position using Integral of Velocity without noise
Pos_noise = array_integral(Vel_measure, delta_t);                           %Position using Integral of Velocity with noise
[Pos_est, Vel_est] = Kalman_filter(sim_time, Vel_measure, delta_t);         %Get Position and Velocity Estimates

Vel_measure1 = Vel_pred + 8*sin(10.1*(2*pi*sim_time));                      %Create Velocity Measure with noise of Sine wave
[Pos_est1, Vel_est1] = Kalman_filter(sim_time, Vel_measure1, delta_t);      %Get Position and Velocity Estimates

[Pos_est2, Vel_est2] = Kalman_filter1(sim_time, Vel_measure, delta_t);      %Kalman Filter with different Q

%% Draw Graph
figure('units', 'pixels', 'pos', [100 100 800 500], 'color', [1, 1, 1], 'Name', 'Random Noise');
    subplot(2,1,1);
    plot(sim_time, Vel_measure, 'k', sim_time, Vel_est, 'r');
    plot_set('측정 / 추정 - 속도', 'Time(s)', 'rad/sec', 15, ...
            [0.0, 0.5, 5, 40, 20, 120]); %Setup Graph
    legend('Vel Measure', 'Vel Est')
   
    subplot(2,1,2);
    plot(sim_time, Pos_ori, 'r', sim_time, Pos_noise, 'b', sim_time, Pos_est, 'k', 'LineWidth',2);
    plot_set('측정 / 노이즈 / 칼만 - 위치', 'Time(s)', 'rad', 15, ...
            [0.0, 0.5, 5, 0, 100, 400]); %Setup Graph
    legend('Pos Ori', 'Pos Noise', 'Pos Est');
    
figure('units', 'pixels', 'pos', [400 100 800 500], 'color', [1, 1, 1], 'Name', 'Sine Noise');
    subplot(2,1,1);
    plot(sim_time, Vel_measure1, 'k', sim_time, Vel_est1, 'r');
    plot_set('측정 / 추정 - 속도', 'Time(s)', 'rad/sec', 15, ...
            [0.0, 0.5, 5, 60, 20, 100]); %Setup Graph
    legend('Vel Measure', 'Vel Est')
   
    subplot(2,1,2);
    plot(sim_time, Pos_ori, 'r', sim_time, Pos_noise, 'b', sim_time, Pos_est1, 'k', 'LineWidth',2);
    plot_set('측정 / 노이즈 / 칼만 - 위치', 'Time(s)', 'rad', 15, ...
            [0.0, 0.5, 5, 0, 100, 400]); %Setup Graph
    legend('Pos Ori', 'Pos Noise', 'Pos Est');


figure('units', 'pixels', 'pos', [900 100 800 500], 'color', [1, 1, 1], 'Name', 'Noise Removal');
    subplot(2,1,1);
    plot(sim_time, Vel_measure, 'k', sim_time, Vel_est2, 'r');
    plot_set('측정 / 추정 - 속도', 'Time(s)', 'rad/sec', 15, ...
            [0.0, 0.5, 5, 40, 20, 120]); %Setup Graph
    legend('Vel Measure', 'Vel Est')
   
    subplot(2,1,2);
    plot(sim_time, Pos_ori, 'r', sim_time, Pos_noise, 'b', sim_time, Pos_est2, 'k', 'LineWidth',2);
    plot_set('측정 / 노이즈 / 칼만 - 위치', 'Time(s)', 'rad', 15, ...
            [0.0, 0.5, 5, 0, 100, 400]); %Setup Graph
    legend('Pos Ori', 'Pos Noise', 'Pos Est');

%% Functions
function p = array_integral(v, dt) %% Integral with return data of Array
    n = length(v);
    p = zeros(1, n);

    for i = 2:n
        f = @(x) interp1(1:i, v(1:i), x, 'linear', 'extrap');
        p(i) = integral(f, 1, i) * dt;  % integral from 1 to i
    end
end

function [Pos_est, Vel_est] = Kalman_filter(sim_time, measure, delta_t) %Kalman Filter
    A = [1 delta_t; 0 1];
    H = [0 1]; %Measuring Speed
    Q = [1 0; 0 3];
    R = 10;
    x = [0; 80];
    P = 5 * eye(2);

    for i = 1:length(sim_time)
        x_hat=A*x;                              % Estimate output
        P_hat=A*P*A'+Q;                         % Estimate Error Covariance
        
        K = P_hat*H'*inv(H*P_hat*H'+R);         % Calculate Kalman Gain
        x = x_hat + K*(measure(i) - H*x_hat);   % Update Output
        P = P_hat - K*H*P_hat;                  % Calculate Error Covariance
    
        Pos_est(i) = x(1);
        Vel_est(i) = x(2);
        P_gain(i) = P(1);
        K_gain(i) = K(1);
    end
end

function [Pos_est, Vel_est] = Kalman_filter1(sim_time, measure, delta_t) %Kalman Filter altered Q
    A = [1 delta_t; 0 1];
    H = [0 1]; %Measuring Speed
    Q = [1 0; 0 0.0001]; %[1 0; 0 3] --> [1 0; 0 0.0001]
    R = 10;
    x = [0; 80];
    P = 5 * eye(2);

    for i = 1:length(sim_time)
        x_hat=A*x;                              % Estimate output
        P_hat=A*P*A'+Q;                         % Estimate Error Covariance
        
        K = P_hat*H'*inv(H*P_hat*H'+R);         % Calculate Kalman Gain
        x = x_hat + K*(measure(i) - H*x_hat);   % Update Output
        P = P_hat - K*H*P_hat;                  % Calculate Error Covariance
    
        Pos_est(i) = x(1);
        Vel_est(i) = x(2);
        P_gain(i) = P(1);
        K_gain(i) = K(1);
    end
end

%Setup graph
function plot_set(Title, XLabel, YLabel, fontsize, tick_param)
    X_min = tick_param(1); X_tick = tick_param(2); X_max = tick_param(3);
    Y_min = tick_param(4); Y_tick = tick_param(5); Y_max = tick_param(6);

    grid on
    axis([X_min X_max Y_min Y_max])
    set(gca, 'XTick', [X_min:X_tick:X_max])
    set(gca, 'YTick', [Y_min:Y_tick:Y_max])

    xlabel(XLabel, 'FontSize', fontsize);
    ylabel(YLabel, 'FontSize', fontsize);
    title(Title, 'FontSize', fontsize);
end