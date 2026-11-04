%为了验证信道稀疏度估计的误差，通过理论计算和实际的仿真对比来观察

%% system parameter
N = 255;
PN_len = 255;
Kmax = 20;
SNR = linspace(0, 30, 15);
test_num = 1000;
alpha = 3.0;
M = 255;
Taps = 20;

%% PN
load pn_all.mat pn420;
PN255 = transpose(pn420(1,420-PN_len+1:420));
PN255 = sqrt(PN_len)*ifft(PN255);
PN255_matrix = zeros(PN_len, PN_len);
for i = 1:PN_len
    PN255_matrix(i:PN_len, i) = PN255(1:PN_len-i+1);
end
Phi = zeros(2*PN_len, PN_len);
for i = 1:PN_len
    Phi(i+PN_len:2*PN_len, i) = PN255(1:PN_len-i+1);
    Phi(i:i+PN_len-1,i) = PN255(1:N);
end

%% K=u[1 taps] simulation of hard threshold method
K_Delta_hard = zeros(test_num, length(SNR));

for i_snr = 1:length(SNR)
    display(i_snr);
    for i_test = 1: test_num
        K = randi([1 Taps]);
        T = randperm(M);
        T = T(1:K);
        X = zeros(N, 1);
        X(T) = random('Normal', 0, 1, K, 1);
        
        Y = Phi*X;
        Y = awgn(Y, SNR(i_snr), 'measured');
        Y = Y(PN_len+1:2*PN_len);
        X_est = ifft(fft(Y, N)./fft(PN255, N), N);
        
        Th = 0.4;
        X_est_filter = (abs(X_est)>Th).*X_est;
        K_Delta_hard(i_test, i_snr) = sum((abs(X_est_filter)>0)~=(abs(X)>0));        
    end
end
K_Delta_hard = mean(K_Delta_hard);
K_Delta_hard = (Taps/2+0.5-K_Delta_hard)/(Taps/2+0.5);


%% K=u[1 taps] simulation of rough method
K_Delta_old = zeros(test_num, length(SNR));

for i_snr = 1:length(SNR)
    display(i_snr);
    for i_test = 1: test_num
        K = randi([1 Taps]);
        T = randperm(M);
        T = T(1:K);
        X = zeros(N, 1);
        X(T) = random('Normal', 0, 1, K, 1);
        
        Y = Phi*X;
        Y = awgn(Y, SNR(i_snr), 'measured');
        Y = Y(PN_len+1:2*PN_len);
        X_est = ifft(fft(Y, N)./fft(PN255, N), N);
        
        Th = 1.9*mean(abs(X_est));
        X_est_filter = (abs(X_est)>Th).*X_est;
        K_Delta_old(i_test, i_snr) = sum((abs(X_est_filter)>0)~=(abs(X)>0));        
    end
end
K_Delta_old = mean(K_Delta_old);
K_Delta_old = (Taps/2+0.5-K_Delta_old)/(Taps/2+0.5);

%% K=u[1 taps] simulation of new method
K_Delta_new = zeros(test_num, length(SNR));
for i_snr = 1:length(SNR)
    display(i_snr);
    for i_test = 1: test_num
        K = randi([1 Taps]);
        T = randperm(M);
        T = T(1:K);
        X = zeros(N, 1);
        X(T) = random('Normal', 0, 1, K, 1);
        
        Y = Phi*X;
        Y = awgn(Y, SNR(i_snr), 'measured');
        Y = Y(PN_len+1:2*PN_len);
        X_est = ifft(fft(Y, N)./fft(PN255, N), N);
        
        Noise = X_est - X;
        Noise_var = var(Noise);
        Noise_mean = mean(Noise);
        Th = alpha*sqrt(Noise_var);
        X_est_filter = (abs(X_est)>Th).*X_est;
        K_Delta_new(i_test, i_snr) = sum((abs(X_est_filter)>0)~=(abs(X)>0));        
    end
end
K_Delta_new = mean(K_Delta_new);
K_Delta_new = (Taps/2+0.5-K_Delta_new)/(Taps/2+0.5);

%% theory K_delta cal

K_Delta_ideal = zeros(1, length(SNR));
for i_snr = 1:length(SNR)
    Sum = 0;
    for K_test = 1:Taps
        Sum = Sum + K_test - K_test*2*(1-normcdf(sqrt(K_test)*3/sqrt(N*(Taps/2)^(0.1*SNR(i_snr))),0,1));
    end    
    K_Delta_ideal(i_snr) = Sum/Taps;
end
display(K_Delta_ideal);
K_Delta_ideal = (Taps/2+0.5-K_Delta_ideal)/(Taps/2+0.5);


%% Plot
figure;
%plot(SNR, mean(K_Delta_10), 'r<-', 'LineWidth', 1.5);
%hold on;
%plot(SNR, K_Delta_10_ideal, 'r<--', 'LineWidth', 1.5);
%hold on;
plot(SNR, K_Delta_hard, 'c+-', 'LineWidth', 1.5);
hold on;
plot(SNR, K_Delta_old, 'b*-', 'LineWidth', 1.5);
hold on;
plot(SNR, K_Delta_ideal, 'ro-', 'LineWidth', 1.5);
hold on;
plot(SNR, K_Delta_new, 'g>--', 'LineWidth', 1.5);
hold on;
grid on;
axis([0 30 0.2 1]);
legend('Hard Threshold', 'Rough Criterions', 'Proposed Pth Threshold', 'Theoretically boundary');
xlabel('SNR (dB)');
ylabel('Detection Probability');
