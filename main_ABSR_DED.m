%% 2014 Adaptive bistable stochastic resonance aided spectrum sensing
%ABSR-DED算法
clearvars
clc
tic
%% 系统基本参数
f_ref = 0.001;          %参考频率
fs_ref = 5;
alpha = 1e6;            %缩放比α

Pf = 0.1;               %虚警率
M = 5000;               %蒙特卡洛仿真次数
%% 信号基本参数
Am = 0.3;       %幅值
fc = 1e3;       %信号频率
fs = 5e6;       %采样率
N = 5e3;        %采样点数
phi = 0;
T = N / fs;     %采样时间

SNR_seq = -30:2:4;
% SNR_seq = -30;
Pd = zeros(size(SNR_seq));
%%
%时域轴
t = (0 : N-1)' / fs;

%无噪信号
s = Am * sin(2*pi*fc*t + phi);

%% 先用各个SNR下的噪声求出虚警率对应的判决阈值(假定噪声功率估计无偏差)
gamma_th_DED = zeros(size(SNR_seq));
for i=1:length(SNR_seq)
    sigma0 = sqrt(Am^2 / 2 / 10^(SNR_seq(i)/10));
    a = 2*pi*f_ref;
    b = a^2 / (2*sigma0^2);
    h = alpha / fs;
    gamma_th_DED(i) = D_ED_threshold(sigma0,N,Pf,[a,b,h]);
end
%% 信号检测仿真
for i=1:length(SNR_seq)
    SNR = SNR_seq(i);      %信噪比
    H_sum = 0;
    gammai = gamma_th_DED(i);
    parfor j=1:M
        %% 产生带噪信号

        %噪声
        sigma = sqrt(Am^2 / 2 / 10^(SNR/10));
        noise = sigma * randn(size(s));

        %带噪信号(接收信号)
        r = s + noise;

        %% 噪声功率估计
        % delta_r = r(1+L : N) - r(1 : N-L);
        % sigma2_MLE = sum(delta_r.^2) / 2 / (N-L); %噪声功率σ2的MLE

        %% 系统最佳参数
        a = 2*pi*f_ref;
        b = a^2 / (2 * sigma^2);
        h = alpha / fs;
        y = Runge_Kutta(a,b,h,r);

        %% D-ED
        %计算信号的检验统计量
        y_mean = mean(y);
        T_y = sum((y-y_mean).^2);

        if T_y > gammai
            H_sum = H_sum + 1;
        end

    end
    Pd(i) = H_sum / M;
end
%% 作图
figure()
plot(SNR_seq,Pd,'r-^','LineWidth',1);
grid on
xlabel('SNR/dB');ylabel('Pd');
axis([-inf inf 0 1]);
legend('ABSR(D-ED)');


toc