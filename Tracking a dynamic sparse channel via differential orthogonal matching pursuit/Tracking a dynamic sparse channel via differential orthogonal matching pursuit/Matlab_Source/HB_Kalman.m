clc;

%% System Parameter
T = 150;
N = 150;
M = 100;
P1 = 0.05;
K = P1*N;
P10 = 0.01;
P01 = K*P10/(N-K);
sigma_h = 1;
sigma_n = 0.05;
SNR = 15;
sigma_a = 0.01;
Pth1 = 5*sigma_n/sqrt(M);
Pth2 = 5*sigma_n/sqrt(M);
OMP_Kmax = 15;
D = 2;

%% Generate Channel
s = zeros(N, T);
a = zeros(N, T);
h = zeros(N, T);
w = zeros(N, T);
s(:,1) = rand(1,N)<P1;
a(:,1) = random('Normal', 0, sigma_h, N, 1);
h(:,1) = s(:,1).*a(:,1);
%figure;
%stem(h(:,1));
for i = 2:T
    K_pre = sum(s(:,i-1));
    for j = 1:N
        if (s(j,i-1)==1)
            if (K_pre>5)
                s(j, i) = rand(1,1)>P10;
            else
                s(j, i) = 1;
            end
        else
            if (K_pre<13) 
                s(j, i) = rand(1,1)<P01;
            else
                s(j, i) = 0;
            end
        end
    end
    %if (i==T/2)
    %    s(:,i) = rand(1,N)<P1;
    %end
    w(:, i) = random('Normal', 0, sigma_a, N, 1);
    a(:, i) = a(:, i-1) + w(:, i);
    h(:, i) = s(:, i).*a(:, i);
    %stem(h(:,i));
end

%% Generate Phi
load pn_all.mat pn420;
PN_len = N;
PN255 = transpose(pn420(1,420-PN_len+1:420));
PN255 = sqrt(PN_len)*ifft(PN255);
PN255_matrix = zeros(PN_len, PN_len);
for i = 1:PN_len
    PN255_matrix(i:PN_len, i) = PN255(1:PN_len-i+1);
end
Phi_l = zeros(2*PN_len, PN_len);
for i = 1:PN_len
    Phi_l(i+PN_len:2*PN_len, i) = PN255(1:PN_len-i+1);
    Phi_l(i:i+PN_len-1,i) = PN255(1:PN_len);
end

Phi = Phi_l(N*2-M+1:N*2, :);
%random('Normal', 0, 1, M, N);
for i = 1:N
    Phi(:,i) = Phi(:,i)/norm(Phi(:,i));
end
%% Generate Y
y = Phi*h;
display(size(Phi));
%y = y + random('Normal', 0, sigma_n, M, T);
y = awgn(y, SNR, 'measured');


%% Kalman sparse filter
h_est_pre = zeros(N, T);
h_est = zeros(N, T);
mse_pre = zeros(N, N, T);
mse = zeros(N, N, T);
h_est(:, 1) = 0;
temp = zeros(1,N);
d_set = find(s(:,1)>0);
temp(d_set) = 1;
mse(:, :, 1) = diag(temp)*sigma_h;
Kmse = zeros(N, M, T);
for i = 2:T
    h_est_pre(:,i) = h_est(:,i-1);
    temp = zeros(1,N);
    d_set = find(s(:,i)>0);
    temp(d_set) = 1;
    mse_pre(:,:,i) = mse(:,:,i-1) + diag(temp)*sigma_a;
    Kmse(:,:,i) = mse_pre(:,:,i)*Phi'*(eye(M)*0.001+Phi*mse_pre(:,:,i)*Phi')^(-1);
    h_est(:,i) = h_est_pre(:,i) + Kmse(:,:,i)*(y(:,i)-Phi*h_est_pre(:,i));
    mse(:,:,i) = (eye(N)-Kmse(:,:,i)*Phi)*mse_pre(:,:,i);
end
mse_Kalman = zeros(1, T);
for i = 1:T
    %h_mse(i) = norm(h(T(i,1:K(i)),i)-h_est(T(i,1:K(i)),i))^2/norm(h(T(i,1:K(i)),i))^2;
    mse_Kalman(i) = norm(h(:,i)-h_est(:,i))^2;
end

%% plot
figure;
semilogy(1:T, mse_Kalman(1:T), 'g', 'LineWidth', 2.5);
hold on;
%axis([1 T 0.5e-3 0.5e-0]);
legend('Kalman');
xlabel('Time');
ylabel('MSE');







